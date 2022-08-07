(in-package #:lla)

(defmethod print-object ((object vector-double) stream)
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

(defmethod print-object ((object vector-float) stream)
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

(defmethod vref ((x vector-float) i)
  (vref@v1 x i))

(defmethod vref ((x vector-double) i)
  (vref@v2 x i))

(defmethod (setf vref) (value (x vector-float) i)
  (setf-vref@v1 value x i))

(defmethod (setf vref) (value (x vector-double) i)
  (setf-vref@v2 value x i))
