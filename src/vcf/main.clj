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
    (s/def ::field     (s/+   (s/alt :unquoted   ::unquoted
                                     :quoted     ::quoted)))
    (s/def ::commaf    (s/cat :comma (set ",")
                              :field ::field))
    (s/def ::commafs   (s/*   ::commaf))
    (s/def ::csv       (s/cat :head ::field
                              :tail ::commafs))
    (s/def ::meta      (s/cat :hashes ::hashes
                              :key    ::key
                              :=      (set "=")
                              :left   (set "<")
                              :csv    ::csv
                              :right  (set ">")))
    (s/def ::<meta>        (s/and string? seqit ::meta))
    (s/def ::contig-line   (meta "contig"))
    (s/def ::contig-line   (meta "ALT"))
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
          [version-line & meta-lines] meta-lines
          [header-lines other-lines] (split-with head? lines)]
      (util/make-map version-line meta-lines header-lines other-lines))))

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
        {:keys [version-line meta-lines header-lines other-lines]} sections
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

(defn parse-meta-line
  "Return [KEY CSV-MAP] from ::<meta> LINE."
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

(defn map-metas
  "Digest all the :meta-lines in VCF into a nested map."
  [vcf]
  (letfn [(kind [[k vs]] [k (map second vs)])
          (id [m] (get m "ID"))
          (noid [m] (dissoc m "ID"))
          (digest [[k ms]] [k (zipmap (map id ms) (map noid ms))])]
    (->> vcf
         section
         :meta-lines
         (map parse-meta-line)
         (remove #{["" {"" ""}]})
         (group-by first)
         (map (comp digest kind))
         (into {})
         walk/keywordize-keys)))

(comment
  (parse vcf)
  (section vcf)
  (conform-lines vcf)
  (time (valid? vcf))
  (map-metas vcf)

  "Not handling these ## meta lines yet."
  ["##Biotin(Bgnd)=Biotin(Bgnd)|Staining|456|292"
   "##Biotin(High)=Biotin(High)|Staining|618|4165"
   "##DNP(Bgnd)=DNP(Bgnd)|Staining|333|266"
   "##DNP(High)=DNP(High)|Staining|12573|355"
   "##Extension(A)=Extension(A)|Extension|25349|405"
   "##Extension(C)=Extension(C)|Extension|1442|8646"
   "##Extension(G)=Extension(G)|Extension|2205|8950"
   "##Extension(T)=Extension(T)|Extension|27882|316"
   "##Hyb(High)=Hyb(High)|Hybridization|2227|8113"
   "##Hyb(Low)=Hyb(Low)|Hybridization|2031|1856"
   "##Hyb(Medium)=Hyb(Medium)|Hybridization|517|5068"
   "##NP(A)=NP(A)|Non-Polymorphic|7725|255"
   "##NP(C)=NP(C)|Non-Polymorphic|568|4775"
   "##NP(G)=NP(G)|Non-Polymorphic|563|4047"
   "##NP(T)=NP(T)|Non-Polymorphic|7892|204"
   "##NSB(Bgnd)Blue=NSB(Bgnd)Blue|Non-SpecificBinding|371|153"
   "##NSB(Bgnd)Green=NSB(Bgnd)Green|Non-SpecificBinding|277|149"
   "##NSB(Bgnd)Purple=NSB(Bgnd)Purple|Non-SpecificBinding|371|161"
   "##NSB(Bgnd)Red=NSB(Bgnd)Red|Non-SpecificBinding|360|153"
   "##Restore=Restore|Restoration|228|305"
   "##String(MM)=String(MM)|Stringency|896|237"
   "##String(PM)=String(PM)|Stringency|13564|257"
   "##TargetRemoval=TargetRemoval|TargetRemoval|1156|208"
   "##analysisVersionNumber=1"
   "##arrayType=Broad_GWAS_supplemental_15061359_A1"
   "##autocallDate=05/23/2018 21:15"
   "##autocallGender=F"
   "##autocallVersion=2.0.0.137"
   "##chipWellBarcode=200557070005_R06C01"
   "##clusterFile=Broad_GWAS_supplemental_15061359_A1.egt"
   "##content=Broad_GWAS_supplemental_15061359_A1.1.2.extended.csv"
   "##expectedGender=Female"
   "##extendedIlluminaManifestVersion=1.2"
   "##extendedManifestFile=Broad_GWAS_supplemental_15061359_A1.1.2.extended.csv"
   "##fileDate=Wed May 23 21:21:12 UTC 2018"
   "##fingerprintGender=Unknown"
   "##genomeBuild=HG19"
   "##imagingDate=3/2/2017 4:24:38 PM"
   "##manifestFile=Broad_GWAS_supplemental_15061359_A1.bpm"
   "##p95Green=4329"
   "##p95Red=8992"
   "##picardVersion=07b46e26eb638116226b10df9f3f653b82b8ea95"
   "##reference=/cromwell_root/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
   "##sampleAlias=NA12878"
   "##scannerName=N370"
   "##source=BPM file"
   "##zcallThresholds=thresholds.7.txt"
   "##zcallVersion=1.0.0.0"])
