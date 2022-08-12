#ifndef LILA_VECTOR_H
#define LILA_VECTOR_H

namespace lila {

SMART(RealSingleVector);

class RealSingleVector_O : public core::CxxObject_O {
  LISP_CLASS(lila, lila_pkg, RealSingleVector_O, "REAL-SINGLE-VECTOR", core::CxxObject_O);
  vector<float> _Value;

public:
  //bool fieldsp() const { return true; };
  //void fields(core::Record_sp node);
  static RealSingleVector_sp make(core::Vaslist_sp args);
};

}; // namespace lila

#endif
