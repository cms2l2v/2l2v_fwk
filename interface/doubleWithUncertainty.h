// <author>Cristóvão B da Cruz e Silva</author>
// <email>c.beirao@cern.ch</email>
// <date>2014-11-24</date>
// <summary>Declaration of a simple class that defines a double with an associated uncertainty</summary>

#ifndef _doubleWithUncertainty_h_
#define _doubleWithUncertainty_h_

#include <iostream>
#include <cmath>

class doubleUnc
{
public:
  doubleUnc(double val, double unc):value_(val),uncertainty2_(unc*unc),defaultUncValue(0){};
  doubleUnc(doubleUnc& val):value_(val.value_),uncertainty2_(val.uncertainty2_),defaultUncValue(0){};

  inline double value() {return value_;};
  inline double uncertainty() {return std::sqrt(uncertainty2_);};
  inline double uncertainty2() {return uncertainty2_;};

  inline double setValue(double value) {return value_ = value;};
  inline double setUncertainty(double value) {return uncertainty2_ = value*value;};
  inline double setUncertainty2(double value) {return uncertainty2_ = std::abs(value);};
  inline int    setDefaultUncValue(int value) {defaultUncValue = 0; if(value == 1) defaultUncValue = 1; if(value == 2) defaultUncValue = 2; return defaultUncValue;};

  double defaultUnc(double currentValue);

  doubleUnc& operator= (doubleUnc& val);
  doubleUnc& operator+=(doubleUnc& val);
  doubleUnc& operator-=(doubleUnc& val);
  doubleUnc& operator*=(doubleUnc& val);
  doubleUnc& operator/=(doubleUnc& val);
  doubleUnc  operator+ (doubleUnc& val);
  doubleUnc  operator- (doubleUnc& val);
  doubleUnc  operator* (doubleUnc& val);
  doubleUnc  operator/ (doubleUnc& val);
  doubleUnc& operator= (double& val);
  doubleUnc& operator+=(double& val);
  doubleUnc& operator-=(double& val);
  doubleUnc& operator*=(double& val);
  doubleUnc& operator/=(double& val);
  doubleUnc  operator+ (double& val);
  doubleUnc  operator- (double& val);
  doubleUnc  operator* (double& val);
  doubleUnc  operator/ (double& val);


  friend std::ostream& operator << (std::ostream &o, doubleUnc& val);

protected:
  double value_;
  double uncertainty2_;
  int defaultUncValue;

private:
};

#endif
