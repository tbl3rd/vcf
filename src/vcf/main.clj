(ns vcf.main
  "Hack VCFs."
  (:require [clojure.string :as str]
            [vcf.util :as util])
  (:gen-class))

(def vcf
  "From gs://broad-gotc-test-storage/annotation_filtration/"
  "./inputs/hg19/negative/200557070005_R06C01.vcf.gz")

(def vcf-version
  "Support only these versions of VCF file."
  #{"VCFv4.2" "VCFv4.3"})

(def INFO-Number
  "1.4.2 Information field format"
  {"A" "a value per alternate allele"
   "R" "a value per possible allele including the reference"
   "G" "a value per posssible genotype"
   "." "number of values varies, is unknown or unbounded"})

(def reserved-INFO-fields
  "from Table 1"
  (do
    [:Field   :Number :Type     :Description]
    [["AA"        "1" "String"  "Ancestral allele"]
     ["AC"        "A" "Integer" "Allele count in genotypes for each ALT in listed order"]
     ["AD"        "R" "Integer" "Total read depth for each allele"]
     ["ADF"       "R" "Integer" "Read depth for each allele on the forward strand"]
     ["ADR"       "R" "Integer" "Read depth for each allele on the reverse strand"]
     ["AF"        "A" "Float"   "Allele frequency for each ALT allele in listed order"]
     ["AN"        "1" "Integer" "Total number of alleles in called genotypes"]
     ["BQ"        "1" "Float"   "RMS base quality"]
     ["CIGAR"     "A" "String"  "How to align an alternate allele to the reference"]
     ["DB"        "0" "Flag"    "dbSNP membership"]
     ["DP"        "1" "Integer" "Combined depth across samples"]
     ["END"       "1" "Integer" "End position for symbolic alleles"]
     ["H2"        "0" "Flag"    "HapMap2 membership"]
     ["H3"        "0" "Flag"    "HapMap3 membership"]
     ["MQ"        "1" "."       "RMS mapping quality"]
     ["MQ0"       "1" "Integer" "Number of MAPQ == 0 reads"]
     ["NS"        "1" "Integer" "Number of samples with data"]
     ["SB"        "4" "Integer" "Strand bias"]
     ["SOMATIC"   "0" "Flag"    "Somatic mutation for cancer genomics"]
     ["VALIDATED" "0" "Flag"    "Validated by follow-up experiment"]
     ["1000G"     "0" "Flag"    "1000 Genomes membership"]]))

(def reserved-genotype-fields
  "1.6.2 Genotype fields"
  (do
    [:Field :Number :Type     :Description]
    [["AD"      "R" "Integer" "Total read depth for each allele"]
     ["ADF"     "R" "Integer" "Read depth for each allele on the forward strand"]
     ["ADR"     "R" "Integer" "Read depth for each allele on the reverse strand"]
     ["DP"      "1" "Integer" "Combined depth across samples"]
     ["EC"      "A" "Integer" "Expected alternate allele counts"]
     ["FT"      "1" "String"  "Filter indicating if this genotype was called"]
     ["GL"      "G" "Float"   "Genotype likelihoods"]
     ["GP"      "G" "Float"   "Genotype posterior probabilities"]
     ["GQ"      "1" "Integer" "Conditional genotype quality"]
     ["GT"      "1" "String"  "Genotype"]
     ["HQ"      "2" "Integer" "Haplotype quality"]
     ["MQ"      "1" "Integer" "RMS mapping quality"]
     ["PL"      "G" "Integer" "Phred-scaled genotype likelihoods rounded to integer"]
     ["PQ"      "1" "Integer" "Phasing quality"]
     ["PS"      "1" "Integer" "Phase set"]]))

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
    (if matches? [:<_> k csv]
        (let [[matches? k csv] (re-matches #"^##([^=]+)=(.*)$" line)]
          (when matches? [:___ k csv])))))

(defn parse-meta
  "Return nil or LINE parsed into a [ID {CSV}] or [ID VALUE]."
  [line]
  (let [string (partial apply str)
        [flag-ignored id csv] (parse-hash-hash line)]
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

;; Add to parsers as they are discovered.
;;
(defn add-a-f-or-i-parser
  "Add to M a parser for FORMAT and INFO meta values."
  [{:strs [Number Type] :as m}]
  (letfn [(->Integer1 [x] (Integer/valueOf x))
          (->Float1   [x] (Float/valueOf   x))]
    (let [parser-for {["1" "Integer"] ->Integer1
                      ["1" "Float"]   ->Float1
                      ["1" "String"]  identity
                      ["A" "Integer"] ->Integer1
                      ["A" "Float"]   ->Float1}]
      (assoc m :parser (parser-for [Number Type] identity)))))

(defn add-parsers
  "Augment FORMAT or INFO metas with an appropriate parser."
  [[meta-key meta-map]]
  [meta-key (if (#{"FORMAT" "INFO"} meta-key)
              (letfn [(add [[k m]] [k (add-a-f-or-i-parser m)])]
                (into {} (map add meta-map)))
              meta-map)])

(defn map-metas
  "Digest all the :meta-lines in VCF into a nested map."
  [vcf]
  (letfn [(fail [s] (throw (new RuntimeException s)))
          (kind [[k vs]] [k (map second vs)])
          (id [m] (get m "ID"))
          (digest [[k vs]]
            [k (cond (every? map? vs)    (zipmap (map id vs) vs)
                     (every? string? vs) (first vs)
                     :else (fail (str "WTF? " (pr-str [k vs]))))])
          (contig-id [line]
            (second (re-find #"^##contig=<ID=([^,]*)," line)))]
    (let [metas (-> vcf section :meta-lines)
          mapped (->> metas
                      (map parse-meta)
                      (group-by first)
                      (map (comp add-parsers digest kind))
                      (into {}))]
      (assoc mapped :contigs (keep contig-id metas)))))

(defn make-variant-mapper
  "A function of LINE using TSV and COLUMNS to return a variant map."
  [tsv columns]
  (let [semi     (util/split-on #";")
        colon    (util/split-on #":")
        k=v      (util/split-on #"=" 2)
        samples  (filter string? columns)
        defaults (zipmap (filter keyword? columns) (repeat identity))
        fixed    (assoc defaults :FILTER semi :FORMAT colon
                        :INFO (fn [info]
                                (into {} (for [[k v] (map k=v (semi info))]
                                           [(keyword k) v]))))
        parser   (apply (partial assoc fixed)
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

(defn count-variants
  "A sorted table of pairs counting variants per contig in VCF."
  [vcf]
  (let [{:keys [contigs]} (map-metas vcf)
        counts (frequencies (map :CHROM (map-variants vcf)))
        chroms (filter (set (keys counts)) contigs)
        table (mapv (fn [chrom] [chrom (counts chrom)]) chroms)]
    (conj table [:TOTAL (apply + (vals counts))])))

(comment
  (map-metas vcf)
  (map-variants vcf)
  (time (count-variants vcf))
  (sort (keys (get (map-metas vcf) "contig")))
  ((get (map-metas vcf) "contig") "MT")
  ((get (map-metas vcf) "contig") "NC_007605")
  ((get (map-metas vcf) "contig") "X")
  ((get (map-metas vcf) "contig") "Y")
  (vals (get (map-metas vcf) "INFO"))
  (distinct (map (fn [m] (select-keys m ["Number" "Type"]))
                 (vals (get (map-metas vcf) "INFO"))))
  (map #(get % "Number") (vals (get (map-metas vcf) "FORMAT")))
  )
