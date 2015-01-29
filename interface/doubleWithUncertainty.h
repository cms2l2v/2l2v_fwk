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
  doubleUnc():value_(0.0),uncertainty2_(0.0){};
  doubleUnc(double val):value_(val),uncertainty2_(0.0){};
  doubleUnc(double val, double unc):value_(val),uncertainty2_(unc*unc),defaultUncValue(0){};
  doubleUnc(const doubleUnc& val):value_(val.value_),uncertainty2_(val.uncertainty2_),defaultUncValue(0){};

//  doubleUnc(const doubleUnc&) = delete;
//  doubleUnc& operator=(const doubleUnc&) = delete;

  inline double value() {return value_;};
  inline double uncertainty() {return std::sqrt(uncertainty2_);};
  inline double uncertainty2() {return uncertainty2_;};

  inline double setValue(double value) {return value_ = value;};
  inline double setUncertainty(double value) {return uncertainty2_ = value*value;};
  inline double setUncertainty2(double value) {return uncertainty2_ = std::abs(value);};
  inline int    setDefaultUncValue(int value) {defaultUncValue = 0; if(value == 1) defaultUncValue = 1; if(value == 2) defaultUncValue = 2; return defaultUncValue;};

  double defaultUnc(double currentValue) const;

  doubleUnc& operator= (const doubleUnc& val);
  doubleUnc& operator+=(const doubleUnc& val);
  doubleUnc& operator-=(const doubleUnc& val);
  doubleUnc& operator*=(const doubleUnc& val);
  doubleUnc& operator/=(const doubleUnc& val);
  doubleUnc  operator+ (const doubleUnc& val) const;
  doubleUnc  operator- (const doubleUnc& val) const;
  doubleUnc  operator* (const doubleUnc& val) const;
  doubleUnc  operator/ (const doubleUnc& val) const;
  doubleUnc& operator= (const double& val);
  doubleUnc& operator+=(const double& val);
  doubleUnc& operator-=(const double& val);
  doubleUnc& operator*=(const double& val);
  doubleUnc& operator/=(const double& val);
  doubleUnc  operator+ (const double& val) const;
  doubleUnc  operator- (const double& val) const;
  doubleUnc  operator* (const double& val) const;
  doubleUnc  operator/ (const double& val) const;


  friend std::ostream& operator << (std::ostream &o, doubleUnc& val);

protected:
  double value_;
  double uncertainty2_;
  int defaultUncValue;

private:
};

#endif
