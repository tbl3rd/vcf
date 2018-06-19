(ns vcf.util
  "Define various useful things."
  (:require [clojure.java.io :as io]
            [clojure.pprint :refer [pprint]])
  (:import [java.util.zip GZIPInputStream]))

(defmacro do-or-nil
  "Value of BODY or nil if it throws."
  [& body]
  `(try (do ~@body)
        (catch Exception x#
          (println x#))))

(defmacro dump
  "Dump [EXPRESSION VALUE] where VALUE is EXPRESSION's value."
  [expression]
  `(let [x# ~expression]
     (do
       (pprint ['~expression x#])
       x#)))

(defmacro make-map
  "Map SYMBOLS as keywords to their values in the environment."
  [& symbols]
  (let [symbols# symbols]
    (zipmap (map keyword symbols) symbols)))

(defn inflate-lines
  "Return a lazy sequence of lines from the gzipped FILE."
  [file]
  (-> file
      io/input-stream
      GZIPInputStream.
      io/reader
      line-seq))
