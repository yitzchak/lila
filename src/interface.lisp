(in-package #:lila)

(defun get-element-tag (args &optional (tag #b00))
  (loop for arg in args
        if (listp arg)
          do (setf tag (get-element-tag arg tag))
        else if (complexp arg)
          do (setf tag (logior #b10 tag))
             (when (or (typep (realpart arg) 'double-float)
                       (typep (imagpart arg) 'double-float))
               (setf tag (logior #b01 tag)))
        else if (typep arg 'double-float)
          do (setf tag (logior #b01 tag)))
  tag)
