(k:includes #~"include/")

(k:systems :lla)

(k:sources :iclasp
           ;#~"src/vector.cc"
           #~"src/lla.cc")

(k:library "cblas" :required t)
