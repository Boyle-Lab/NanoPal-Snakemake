(ql:quickload '(:losh :iterate :jarl :conserve))

;; Script for taking the brain experiments sheet TSV and generating the list of
;; samples on Turbo to pass to the pipeline.  Might need something similar in
;; the future so I'm saving it here.

(defpackage :samples
  (:use :cl :losh :iterate))

(in-package :samples)

(defun target-id (key)
  (substitute #\- #\space (format nil "~A_~A_~A" (first key) (second key) (third key))))

(defun target-path (entry)
  (let ((pair (drop 3 entry)))
    (unless (string= "" (first pair))
      (format nil "~A/~A.sorted.bam" (first pair) (second pair)))))

(defparameter *data*
  (with-open-file (f "samples.tsv")
    (let ((conserve:*delimiter* #\tab))
      (_ f
        conserve:read-rows
        (remove-if (lambda (row) (string= "" (fourth row))) _)
        (group-by (curry #'take 3) _ :test 'equal)
        (mutate-hash-values (curry #'mapcar #'target-path) _)))))

(defparameter *individual*
  (_ *data*
    (iterate (for (k v) :in-hashtable _)
             (collect-hash ((target-id k) v) :test 'equal))
    (mutate-hash-values (rcurry #'coerce 'vector) _)))

(defparameter *bulk*
  (_ *data*
    alexandria:hash-table-alist
    (group-by #'caar _ :test 'equal)
    (mutate-hash-values (curry #'mapcar #'cdr) _)
    (mutate-hash-values (curry #'reduce #'nconc) _)
    (iterate (for (k v) :in-hashtable _)
             (collect-hash ((format nil "bulk_~A" k) v) :test 'equal))
    (mutate-hash-values (rcurry #'coerce 'vector) _)))



;; (pbcopy (jarl:print *individual* nil))

;; (pbcopy (jarl:print *bulk* nil))
