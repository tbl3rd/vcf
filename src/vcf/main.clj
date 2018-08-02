(ns vcf.main
  "Hack VCFs."
  (:require [clojure.string :as str]
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

(defn make-variant-mapper
  "A function of LINE using TSV and COLUMNS to return a variant map."
  [tsv columns]
  (let [semi  (util/split-on #";")
        colon (util/split-on #":")
        k=v   (util/split-on #"=" 2)
        samples (filter string? columns)
        defaults (zipmap (filter keyword? columns) (repeat identity))
        fixed (assoc
               defaults :FILTER semi :FORMAT colon
               :INFO (fn [info]
                       (into {} (for [[k v] (map k=v (semi info))]
                                  [(keyword k) v]))))
        parser (apply (partial assoc fixed)
                      (interleave samples (repeat colon)))]
    (letfn [(parse [m [k v]] (assoc m k ((parser k) v)))]
      (fn [line]
        (let [m (reduce parse {} (zipmap columns (tsv line)))
              formats (partial zipmap (map keyword (:FORMAT m)))]
          (map (fn [sample] (-> m (update sample formats) (dissoc :FORMAT)))
               samples))))))

(defn map-variants
  "Map variant lines from VCF with headers distributed over samples."
  [vcf]
  (let [tsv (util/split-on #"\t")
        {:keys [header-lines other-lines]} (section vcf)
        headers (-> header-lines first (subs 1) tsv)
        [keywords samples] (split-at 9 headers)
        columns (concat (map keyword keywords) samples)]
    (mapcat (make-variant-mapper tsv columns) other-lines)))

(comment
  (section vcf)
  (map-metas vcf)
  (map-variants vcf)
  (map #(get % "Number") (vals (get (map-metas vcf) "INFO")))
  (map #(get % "Number") (vals (get (map-metas vcf) "FORMAT")))
  (frequencies (map :CHROM (map-variants vcf)))
  )
