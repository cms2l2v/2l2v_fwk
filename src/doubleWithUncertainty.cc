// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-11-24</date>
// <summary>Definition of a simple class that defines a double with an associated uncertainty</summary>

#include "UserCode/llvv_fwk/interface/doubleWithUncertainty.h"

std::ostream& operator << (std::ostream &o, doubleUnc& val)
{
  return o << val.value_ << " +- " << std::sqrt(val.uncertainty2_);
}

doubleUnc& doubleUnc::operator= (doubleUnc& val)
{
  value_ = val.value_;
  uncertainty2_ = val.uncertainty2_;

  return *this;
}

doubleUnc& doubleUnc::operator+=(doubleUnc& val)
{
  value_ += val.value_;
  uncertainty2_ += val.uncertainty2_;

  return *this;
}

doubleUnc& doubleUnc::operator-=(doubleUnc& val)
{
  value_ -= val.value_;
  uncertainty2_ += val.uncertainty2_;

  return *this;
}

doubleUnc& doubleUnc::operator*=(doubleUnc& val)
{
  uncertainty2_ = val.value_*val.value_*uncertainty2_ + value_*value_*val.uncertainty2_;
  value_ *= val.value_;

  return *this;
}

doubleUnc& doubleUnc::operator/=(doubleUnc& val)
{
  uncertainty2_ = uncertainty2_/(val.value_*val.value_) + (val.uncertainty2_/(val.value_*val.value_)) * ((value_*value_)/(val.value_*val.value_));
  value_ /= val.value_;

  return *this;
}

doubleUnc  doubleUnc::operator+ (doubleUnc& val)
{
  doubleUnc retVal(*this);

  return retVal += val;
}

doubleUnc  doubleUnc::operator- (doubleUnc& val)
{
  doubleUnc retVal(*this);

  return retVal -= val;
}

doubleUnc  doubleUnc::operator* (doubleUnc& val)
{
  doubleUnc retVal(*this);

  return retVal *= val;
}

doubleUnc  doubleUnc::operator/ (doubleUnc& val)
{
  doubleUnc retVal(*this);

  return retVal /= val;
}

doubleUnc& doubleUnc::operator= (double& val)
{
  value_ = val;
  uncertainty2_ = defaultUnc(val);

  return *this;
}

doubleUnc& doubleUnc::operator+=(double& val)
{
  doubleUnc tempVal(val, defaultUnc(val));
  return (*this) += tempVal;
}

doubleUnc& doubleUnc::operator-=(double& val)
{
  doubleUnc tempVal(val, defaultUnc(val));
  return (*this) -= tempVal;
}

doubleUnc& doubleUnc::operator*=(double& val)
{
  doubleUnc tempVal(val, defaultUnc(val));
  return (*this) *= tempVal;
}

doubleUnc& doubleUnc::operator/=(double& val)
{
  doubleUnc tempVal(val, defaultUnc(val));
  return (*this) /= tempVal;
}

doubleUnc  doubleUnc::operator+ (double& val)
{
  doubleUnc retVal(val, defaultUnc(val));

  return retVal += *this;
}

doubleUnc  doubleUnc::operator- (double& val)
{
  doubleUnc tempVal(val, defaultUnc(val));
  doubleUnc retVal = *this;
  return retVal -= tempVal;
}

doubleUnc  doubleUnc::operator* (double& val)
{
  doubleUnc retVal(val, defaultUnc(val));

  return retVal *= *this;
}

doubleUnc  doubleUnc::operator/ (double& val)
{
  doubleUnc tempVal(val, defaultUnc(val));
  doubleUnc retVal = *this;

  return retVal /= tempVal;
}

double doubleUnc::defaultUnc(double currentValue)
{
  switch(defaultUncValue)
  {
  case 0:
    return 0;
    break;
  case 1:
    return 1;
    break;
  case 2:
    return std::sqrt(currentValue);
    break;
  }

  return 0;
}
