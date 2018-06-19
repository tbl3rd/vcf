#!/usr/bin/env boot

(set-env!
 :resource-paths #{"src"}
 :target-path "target"
 :repositories
 '[["clojars"       {:url "https://repo.clojars.org/"}]
   ["maven-central" {:url "http://repo.maven.apache.org/maven2/"}]]
 :dependencies '[[org.clojure/clojure    "1.9.0"]
                 [org.clojure/spec.alpha "0.2.168"]
                 [org.clojure/test.check "0.9.0" :scope "test"]])

;; So boot.lein can pick up the project name and version.
;;
(task-options! pom {:project 'vcf
                    :version "0"
                    :developers {"Tom Lyons" "tbl@broadinstitute.org"}})

(deftask link
  "Link ./vcf to ./build.boot to run as a script."
  []
  (with-pre-wrap fileset
    (when-not (.exists (clojure.java.io/file "." "vcf"))
      (dosh "ln" "-s" "./build.boot" "vcf"))
    fileset))

(deftask build
  "Build this."
  []
  (comp (link)
        (aot :namespace '#{vcf.main})
        (uber)
        (pom)
        (jar :main 'vcf.main
             :manifest {"Description" "Hack VCFs."})
        (target)))
