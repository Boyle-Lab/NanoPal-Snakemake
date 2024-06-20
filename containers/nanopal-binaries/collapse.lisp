(ql:quickload '(:cl-interval :iterate :losh :cl-ppcre :str :adopt))

(defpackage :collapse
  (:use :cl :iterate :losh))

(in-package :collapse)

(defun run (path)
  (iterate
    (with db = (make-hash-table :test 'equal))
    (for line :in-file path :using #'read-line)
    (for row = (str:words line))
    (for id = (fourth row))
    (push row (gethash id db))
    (finally (iterate
               (for (nil rows) :in-hashtable db)
               (for kinds = (nreverse (mapcar #'fifth rows)))
               (for (chr start end id) = (first rows))
               (write-line (str:join #\tab (list chr start end id (str:join #\/ kinds))))))))

(defun toplevel ()
  (run (second (adopt:argv))))

(defun build ()
  (sb-ext:save-lisp-and-die "collapse"
    :toplevel #'toplevel
    :executable t :compression t :save-runtime-options t))

;; (run "containers/nanopal-binaries/bin/input_rm_cluster.txt")
