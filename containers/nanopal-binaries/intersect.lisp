(ql:quickload '(:cl-interval :iterate :losh :cl-ppcre :str :adopt))

(defpackage :intersect
  (:use :cl :iterate :losh))

(in-package :intersect)

(defparameter *re* (ppcre:create-scanner "[ \\t]+"))

(defun parse-line (line)
  (destructuring-bind (chr start end tail)
      (ppcre:split *re* line :limit 4)
    (list chr (- (parse-integer start) 50) (+ (parse-integer end) 50) tail)))

(defun parse-file (path)
  (iterate
    (for line :in-file path :using #'read-line)
    (collect (parse-line line))))

(defstruct single-read
  chr
  start
  end
  tail)

(defstruct (interval (:include interval:interval))
  (members nil))

(defun print-chunk (chr start end sep tail)
  (princ chr) (princ #\tab)
  (princ start) (princ #\tab)
  (princ end) (princ #\tab)
  (princ sep)
  (princ tail))

(defun run (path1 path2)
  (let ((data1 (parse-file path1))
        (data2 (parse-file path2))
        (db (make-hash-table :test 'equal)))
    (loop :for (chr start end tail) :in data2
          :for r = (make-single-read :chr chr :start start :end end :tail tail)
          :for tree = (alexandria:ensure-gethash chr db (interval:make-tree))
          :for interval = (interval:insert tree (make-interval :start start :end end :members nil))
          :do (push r (interval-members interval)))
    (loop :for (chr start end tail) :in data1
          :for tree = (alexandria:ensure-gethash chr db (interval:make-tree))
          :for intersects = (interval:find-all tree (cons start end))
          :do (if intersects
                (dolist (i intersects)
                  (dolist (r (interval-members i))
                    (print-chunk chr start end #\space tail)
                    (princ #\tab)
                    (print-chunk (single-read-chr r)
                                 (single-read-start r)
                                 (single-read-end r)
                                 #\tab
                                 (single-read-tail r))
                    (terpri)))
                (progn (print-chunk chr start end #\space tail)
                       (terpri))))))

(defun toplevel ()
  (destructuring-bind (path1 path2) (rest (adopt:argv))
    (run path1 path2)))

(defun build (dest)
  (sb-ext:save-lisp-and-die dest
    :toplevel #'toplevel
    :executable t :compression t :save-runtime-options t))
