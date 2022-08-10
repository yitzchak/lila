(in-package #:lila)

(shadow '("VECTORP"))

(defgeneric vectorp (x)
  (:method (x)
    (declare (ignore x))
    nil))

(defgeneric dimension (x))

(defgeneric l1-norm (x))

(defgeneric l2-norm (x))

(defgeneric l2-norm-sqr (x))

(defgeneric vref (x i))

(defgeneric (setf vref) (new-value v i))

(defun dot (x y)
  (inner x y))

(export '(dimension dot l1-norm l2-norm l2-norm-sqr vref vectorp))
