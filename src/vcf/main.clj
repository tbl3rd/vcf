(ns vcf.main
  "Hack VCFs."
  (:require [clojure.spec.alpha :as s]
            [clojure.string :as str]
            [clojure.walk :as walk]
            [vcf.util :as util])
  (:gen-class))

(def vcf-version
  "Support only these versions of VCF file."
  #{"VCFv4.2" "VCFv4.3"})

(def separators
  "Map column key to field separator as a regex."
  {"FORMAT" #":"
   "INFO"   #";"})

(def vcf
  "From gs://broad-gotc-test-storage/annotation_filtration/"
  "./inputs/hg19/negative/200557070005_R06C01.vcf.gz")

;; Sketch a loose spec of a VCF file.
;; See: https://samtools.github.io/hts-specs/VCFv4.3.pdf
;;
(let [lower "abcdefghijklmnopqrstuvwxyz"
      upper (str/upper-case lower)
      digit "0123456789"
      other "(_)"
      seqit (s/conformer seq (partial apply str))]
  (letfn [(is? [k v] (fn [x] (v (apply str (k x)))))
          (meta [k] (s/and string? seqit ::meta (is? :key #{k})))]
    (s/def ::dq        (set "\""))
    (s/def ::not-dq    (complement (set "\"")))
    (s/def ::whatever  (s/* (constantly true)))
    (s/def ::hash      (set "#"))
    (s/def ::hashes    (s/cat :hash0 ::hash :hash1 ::hash))
    (s/def ::key       (s/+ (set (concat lower upper digit other))))
    (s/def ::not-hash  (complement (set "#")))
    (s/def ::tab       #{\tab})
    (s/def ::not-tab   (complement #{\tab}))
    (s/def ::key-value (s/cat :key ::key := (set "=") :value ::whatever))
    (s/def ::column    (s/+ ::not-tab))
    (s/def ::header    (s/cat :tab ::tab :column ::column))
    (s/def ::hash-line (s/and string? seqit
                              (s/cat :hashes ::hashes
                                     :key    ::key
                                     :=      (set "=")
                                     :value  ::whatever)))
    (s/def ::form-line (s/and ::hash-line
                              (is? :key   #{"fileformat"})
                              (is? :value vcf-version)))
    (s/def ::assembly  (s/and ::hash-line
                              (is? :key #{"assembly"})))
    (s/def ::head-line (s/and string? seqit
                              (s/cat :hash    ::hash
                                     :chrom   ::column
                                     :headers (s/+ ::header))
                              (is? :chrom #{"CHROM"})))
    (s/def ::body-line (s/and string? seqit
                              (s/cat :not-hash ::not-hash
                                     :whatever ::whatever)))
    (s/def ::vcf-line  (s/or  :form       ::form-line
                              :meta       ::hash-line
                              :head       ::head-line
                              :body       ::body-line))
    (s/def ::vcf-file  (s/cat :form-line  ::form-line
                              :hash-lines (s/* ::hash-line)
                              :head-line  ::head-line
                              :body-lines (s/* (s/or :head ::head-line
                                                     :body ::body-line))))
    (s/def ::quoted    (s/cat :dq    ::dq
                              :value (s/+ ::not-dq)
                              :dq    ::dq))
    (s/def ::unquoted  (s/+   (complement (set ",\""))))
    (s/def ::field     (s/+   (s/alt :unquoted ::unquoted
                                     :quoted   ::quoted)))
    (s/def ::commaf    (s/cat :comma (set ",")
                              :field ::field))
    (s/def ::commafs   (s/*   ::commaf))
    (s/def ::csv       (s/cat :head ::field
                              :tail ::commafs))
    (s/def ::meta-kv   (s/and string? seqit
                              (s/cat  :hashes ::hashes
                                      :key-value ::key-value)))
    (s/def ::meta      (s/cat :hashes ::hashes
                              :key    ::key
                              :=      (set "=")
                              :left   (set "<")
                              :csv    ::csv
                              :right  (set ">")))
    (s/def ::<meta>        (s/and string? seqit ::meta))
    (s/def ::contig-line   (meta "contig"))
    (s/def ::alt-line      (meta "ALT"))
    (s/def ::filter-line   (meta "FILTER"))
    (s/def ::format-line   (meta "FORMAT"))
    (s/def ::info-line     (meta "INFO"))
    (s/def ::meta-line     (meta "META"))
    (s/def ::pedigree-line (meta "PEDIGREE"))
    (s/def ::sample-line   (meta "SAMPLE"))))

(defn valid?
  "True when first 9999 lines of VCF is a valid VCF file.
  This takes about a minute to run now!"
  [vcf]
  (s/valid? ::vcf-file (->> vcf
                            util/inflate-lines
                            (take 9999))))

(defn conform-lines
  "Lazily parse the lines of VCF by applying s/conform."
  [vcf]
  (letfn [(parse [line] (s/conform ::vcf-line line))]
    (->> vcf
         util/inflate-lines
         (map parse))))

(defn section
  "Return a map of the sections of the VCF."
  [vcf]
  (letfn [(meta? [line] (str/starts-with? line "##"))
          (head? [line] (str/starts-with? line "#"))]
    (let [[meta-lines lines] (->> vcf
                                  util/inflate-lines
                                  (split-with meta?))
          [header-lines other-lines] (split-with head? lines)]
      (util/make-map meta-lines header-lines other-lines))))

(defn validate-lines
  "Throw when LINES are not valid for a VCF file."
  [vcf version-line head-line header-lines]
  (letfn [(fail [s] (throw (new RuntimeException s)))
          (check [b s] (when-not b (fail s)))]
    (check (s/valid? ::form-line version-line)
           (str vcf " is not a version '" vcf-version "' VCF file."))
    (check (s/valid? ::head-line (first header-lines))
           (str "No #CHROM header line in " vcf))
    (check (every? #{head-line} header-lines)
           (str "Conflicting header lines in " vcf))))

(defn parse
  "Return a map representing the content of VCF file."
  [vcf]
  (let [sections (section vcf)
        {:keys [meta-lines header-lines other-lines]} sections
        version-line (first meta-lines)
        head-line (first header-lines)]
    (validate-lines vcf version-line head-line header-lines)
    (let [version (->> version-line
                       (s/conform ::form-line)
                       :value
                       (apply str))
          {:keys [chrom headers]} (s/conform ::head-line head-line)
          columns (->> headers
                       (map :column)
                       (into [chrom])
                       (map (partial apply str)))
          variants (->> other-lines
                        (remove (fn [line] (= head-line line)))
                        (map (fn [line] (str/split line #"\t"))))]
      (util/make-map version columns variants))))

(defn parse-meta-line-slow
  "Return [KEY CSV-MAP] from ::<meta> LINE using spec."
  [line]
  (util/dump line)
  (letfn [(fieldit [field] (->> field
                                (s/unform ::field)
                                (s/conform ::key-value)
                                ((juxt :key :value))
                                (mapv (partial apply str))))]
    (let [{:keys [key csv]} (s/conform ::<meta> line)
          {:keys [head tail]} csv
          fields (into [head] (map :field tail))]
      [(apply str key) (into {} (map fieldit fields))])))

(defn map-metas-slow
  "Digest all the :meta-lines in VCF into a nested map."
  [vcf]
  (letfn [(kind [[k vs]] [k (map second vs)])
          (id [m] (get m "ID"))
          (noid [m] (dissoc m "ID"))
          (digest [[k ms]] [k (zipmap (map id ms) (map noid ms))])]
    (->> vcf
         section
         :meta-lines
         (map parse-meta-line-slow)
         (remove #{["" {"" ""}]})
         (group-by first)
         (map (comp digest kind))
         (into {})
         walk/keywordize-keys)))

(defn parse-hash-hash
  "Parse a ##ID=CSV or ##ID=<CSV> LINE returning [ID CSV] or nil."
  [line]
  (let [[matches? k csv] (re-matches #"^##([^<]+)=<(.*)>$" line)]
    (if matches? [k csv]
        (let [[matches? k csv] (re-matches #"^##([^=]+)=(.*)$" line)]
          (when matches? [k csv])))))

(defn parse-meta
  "Return nil or LINE parsed into a [ID {CSV}] or [ID VALUE]."
  [line]
  (let [string (partial apply str)
        [id csv] (parse-hash-hash line)]
    (when id
      (loop [m {} dq nil k [] eq nil v [] in csv]
        (if-let [c (first in)]
          (case c
            \" (if dq
                 (recur m nil k eq v (rest in))
                 (recur m c k eq v (rest in)))
            \, (if dq
                 (recur m dq k  eq  (conj v c) (rest in))
                 (recur (assoc m (string k) (string v)) dq [] nil [] (rest in)))
            \= (if dq
                 (if eq
                   (recur m dq k eq (conj v c) (rest in))
                   (recur m dq (conj k c) eq v (rest in)))
                 (recur m dq k c v (rest in)))
            (if eq
              (recur m dq k eq (conj v c) (rest in))
              (recur m dq (conj k c) eq v (rest in))))
          [id (if (seq v)
                (assoc m (string k) (string v))
                (string k))])))))

(defn map-metas
  "Digest all the :meta-lines in VCF into a nested map."
  [vcf]
  (letfn [(fail [s] (throw (new RuntimeException s)))
          (kind [[k vs]] [k (map second vs)])
          (id [m] (get m "ID"))
          (digest [[k vs]]
            [k (cond
                 (every? map? vs)    (zipmap (map id vs) vs)
                 (every? string? vs) (first vs)
                 :else (fail (str "WTF? " (pr-str [k vs]))))])]
    (->> vcf
         section
         :meta-lines
         (map parse-meta)
         (group-by first)
         (map (comp digest kind))
         (into {}))))

(comment
  (parse vcf)
  (section vcf)
  (conform-lines vcf)
  (time (valid? vcf))
  (map-metas vcf)
  )
