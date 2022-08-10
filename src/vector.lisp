(in-package #:lila)

(defun read-vector (stream char n)
  (declare (ignore char n))
  (let* ((x (read stream t nil t))
         (complexp (some #'complexp x))
         (doublep (some (lambda (y)
                          (or (typep y 'double-float)
                              (and (complexp y)
                                   (or (typep (realpart y) 'double-float)
                                       (typep (imagpart y) 'double-float)))))
                        x)))
    (apply (cond ((and complexp doublep) #'complex-double-vector)
                 (complexp #'complex-single-vector)
                 (doublep #'real-double-vector)
                 (t #'real-single-vector))
           x)))
  
(set-dispatch-macro-character #\# #\V #'read-vector)

(defmethod vectorp ((x real-single-vector))
  (declare (ignore x))
  t)

(defmethod vectorp ((x real-double-vector))
  (declare (ignore x))
  t)

(defmethod vectorp ((x complex-single-vector))
  (declare (ignore x))
  t)

(defmethod vectorp ((x complex-double-vector))
  (declare (ignore x))
  t)

(defmethod dimension ((x real-single-vector))
  (dimension@r1v x))

(defmethod dimension ((x real-double-vector))
  (dimension@r2v x))

(defmethod dimension ((x complex-single-vector))
  (dimension@c1v x))

(defmethod dimension ((x complex-double-vector))
  (dimension@c2v x))

(defmethod print-object ((object real-double-vector) stream)
  (if *print-readably*
      (loop for i below (dimension object)
            finally (write-char #\) stream)
            when (cl:zerop i)
              do (write-string "#V(" stream)
            else
              do (write-char #\Space stream)
            do (write (vref object i) :stream stream))
      (print-unreadable-object (object stream :type t)
        (loop for i below (dimension object)
              unless (cl:zerop i)
                do (write-char #\Space stream)
              do (write (vref object i) :stream stream))))
  object)

(defmethod print-object ((object real-single-vector) stream)
  (if *print-readably*
      (loop for i below (dimension object)
            finally (write-char #\) stream)
            when (cl:zerop i)
              do (write-string "#V(" stream)
            else
              do (write-char #\Space stream)
            do (write (vref object i) :stream stream))
      (print-unreadable-object (object stream :type t)
        (loop for i below (dimension object)
              unless (cl:zerop i)
                do (write-char #\Space stream)
              do (write (vref object i) :stream stream))))
  object)

(defmethod print-object ((object complex-double-vector) stream)
  (if *print-readably*
      (loop for i below (dimension object)
            finally (write-char #\) stream)
            when (cl:zerop i)
              do (write-string "#V(" stream)
            else
              do (write-char #\Space stream)
            do (write (vref object i) :stream stream))
      (print-unreadable-object (object stream :type t)
        (loop for i below (dimension object)
              unless (cl:zerop i)
                do (write-char #\Space stream)
              do (write (vref object i) :stream stream))))
  object)

(defmethod print-object ((object complex-single-vector) stream)
  (if *print-readably*
      (loop for i below (dimension object)
            finally (write-char #\) stream)
            when (cl:zerop i)
              do (write-string "#V(" stream)
            else
              do (write-char #\Space stream)
            do (write (vref object i) :stream stream))
      (print-unreadable-object (object stream :type t)
        (loop for i below (dimension object)
              unless (cl:zerop i)
                do (write-char #\Space stream)
              do (write (vref object i) :stream stream))))
  object)

(defmethod vref ((x real-single-vector) i)
  (vref@r1v x i))

(defmethod vref ((x real-double-vector) i)
  (vref@r2v x i))

(defmethod vref ((x complex-single-vector) i)
  (vref@c1v x i))

(defmethod vref ((x complex-double-vector) i)
  (vref@c2v x i))

(defmethod (setf vref) (value (x real-single-vector) i)
  (setf-vref@r1v value x i))

(defmethod (setf vref) (value (x real-double-vector) i)
  (setf-vref@r2v value x i))

(defmethod (setf vref) (value (x complex-single-vector) i)
  (setf-vref@c1v value x i))

(defmethod (setf vref) (value (x complex-double-vector) i)
  (setf-vref@c2v value x i))

(defmethod l1-norm ((x real-single-vector))
  (l1-norm@r1v x))

(defmethod l1-norm ((x real-double-vector))
  (l1-norm@r2v x))

(defmethod l1-norm ((x complex-single-vector))
  (l1-norm@c1v x))

(defmethod l1-norm ((x complex-double-vector))
  (l1-norm@c2v x))

(defmethod l2-norm ((x real-single-vector))
  (l2-norm@r1v x))

(defmethod l2-norm ((x real-double-vector))
  (l2-norm@r2v x))

(defmethod l2-norm ((x complex-single-vector))
  (l2-norm@c1v x))

(defmethod l2-norm ((x complex-double-vector))
  (l2-norm@c2v x))

(defmethod l2-norm-sqr ((x real-single-vector))
  (l2-norm-sqr@r1v x))

(defmethod l2-norm-sqr ((x real-double-vector))
  (l2-norm-sqr@r2v x))

(defmethod l2-norm-sqr ((x complex-single-vector))
  (l2-norm-sqr@c1v x))

(defmethod l2-norm-sqr ((x complex-double-vector))
  (l2-norm-sqr@c2v x))
