(in-package #:lila)

(defun (setf vref) (new-value v i)
  (setf-vref v new-value i))
