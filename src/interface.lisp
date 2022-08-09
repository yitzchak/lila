(in-package #:lla)

(defgeneric inner (x y)
  (:method ((x vector-float) (y vector-float))
    (inner@v1@v1 x y))
  (:method ((x vector-double) (y vector-double))
    (inner@v2@v2 x y)))
