use std::cmp::{Eq, PartialEq};
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::geometry::dq::DualQuaternion;

// Equality: DualQuaternion
impl PartialEq<DualQuaternion> for DualQuaternion {
  fn eq(&self, other: &DualQuaternion) -> bool {
    let q1 = self;
    let q2 = other;
    q1.p == q2.p && q1.d == q2.d
  }
}

// Equality: RHS scalar
impl PartialEq<f64> for DualQuaternion {
  fn eq(&self, other: &f64) -> bool {
    let q1 = self;
    let q2 = DualQuaternion::new_real_quaternion(*other);
    q1.p == q2.p && q1.d == q2.d
  }
}

// Equality: LHS scalar
impl PartialEq<DualQuaternion> for f64 {
  fn eq(&self, other: &DualQuaternion) -> bool {
    let q1 = DualQuaternion::new_real_quaternion(*self);
    let q2 = other;
    q1.p == q2.p && q1.d == q2.d
  }
}

impl Eq for DualQuaternion {}

// Addition: DualQuaternion
impl Add<DualQuaternion> for DualQuaternion {
  type Output = DualQuaternion;

  fn add(self, q2: DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p + q2.p;
    let d = q1.d + q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

// Addition: DualQuaternion Reference
impl Add<&DualQuaternion> for DualQuaternion {
  type Output = DualQuaternion;

  fn add(self, q2: &DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p + q2.p;
    let d = q1.d + q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

impl Add<DualQuaternion> for &DualQuaternion {
  type Output = DualQuaternion;

  fn add(self, q2: DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p + q2.p;
    let d = q1.d + q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

impl Add<&DualQuaternion> for &DualQuaternion {
  type Output = DualQuaternion;

  fn add(self, q2: &DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p + q2.p;
    let d = q1.d + q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

// Addition: RHS scalar
impl Add<f64> for DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn add(self, rhs: f64) -> DualQuaternion {
    self + DualQuaternion::new_real_quaternion(rhs)
  }
}

// Addition: LHS scalar
impl Add<DualQuaternion> for f64 {
  type Output = DualQuaternion;

  #[inline]
  fn add(self, rhs: DualQuaternion) -> DualQuaternion {
    DualQuaternion::new_real_quaternion(self) + rhs
  }
}

// Subtraction: DualQuaternion
impl Sub<DualQuaternion> for DualQuaternion {
  type Output = DualQuaternion;

  fn sub(self, q2: DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p - q2.p;
    let d = q1.d - q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

// Subtraction: DualQuaternion Reference
impl Sub<&DualQuaternion> for DualQuaternion {
  type Output = DualQuaternion;

  fn sub(self, q2: &DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p - q2.p;
    let d = q1.d - q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

impl Sub<DualQuaternion> for &DualQuaternion {
  type Output = DualQuaternion;

  fn sub(self, q2: DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p - q2.p;
    let d = q1.d - q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

impl Sub<&DualQuaternion> for &DualQuaternion {
  type Output = DualQuaternion;

  fn sub(self, q2: &DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p - q2.p;
    let d = q1.d - q2.d;

    DualQuaternion::from_parts(p, d)
  }
}

// Substraction: RHS scalar
impl Sub<f64> for DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn sub(self, rhs: f64) -> DualQuaternion {
    self - DualQuaternion::new_real_quaternion(rhs)
  }
}

// Substraction: LHS scalar
impl Sub<DualQuaternion> for f64 {
  type Output = DualQuaternion;

  #[inline]
  fn sub(self, rhs: DualQuaternion) -> DualQuaternion {
    DualQuaternion::new_real_quaternion(self) - rhs
  }
}

// Multiplication: DualQuaternion
impl Mul<DualQuaternion> for DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn mul(self, q2: DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p * q2.p;
    let d = q1.p * q2.d + q1.d * q2.p;

    DualQuaternion::from_parts(p, d)
  }
}

// Multiplaction: DualQuaternion Reference
impl Mul<&DualQuaternion> for DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn mul(self, q2: &DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p * q2.p;
    let d = q1.p * q2.d + q1.d * q2.p;

    DualQuaternion::from_parts(p, d)
  }
}

impl Mul<DualQuaternion> for &DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn mul(self, q2: DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p * q2.p;
    let d = q1.p * q2.d + q1.d * q2.p;

    DualQuaternion::from_parts(p, d)
  }
}

impl Mul<&DualQuaternion> for &DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn mul(self, q2: &DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let p = q1.p * q2.p;
    let d = q1.p * q2.d + q1.d * q2.p;

    DualQuaternion::from_parts(p, d)
  }
}

// Multiplaction: RHS scalar
impl Mul<f64> for DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn mul(self, rhs: f64) -> DualQuaternion {
    DualQuaternion::from_parts(self.p * rhs, self.d * rhs)
  }
}

// Multiplaction: LHS scalar
impl Mul<DualQuaternion> for f64 {
  type Output = DualQuaternion;

  #[inline]
  fn mul(self, rhs: DualQuaternion) -> DualQuaternion {
    DualQuaternion::from_parts(self * rhs.p, self * rhs.d)
  }
}

// Division: DualQuaternion
impl Div<DualQuaternion> for DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn div(self, rhs: DualQuaternion) -> DualQuaternion {
    let q1 = &self;
    let q2 = &rhs;

    let a = q1.p;
    let b = q1.d;
    let c = q2.p;
    let d = q2.d;

    DualQuaternion::from_parts(
      a.right_div(&c).unwrap(),
      (b * c - a * d).right_div(&c.powf(2.)).unwrap(),
    )
  }
}

// Negation: DualQuaternion
impl Neg for DualQuaternion {
  type Output = DualQuaternion;

  #[inline]
  fn neg(self) -> DualQuaternion {
    let p = -self.p;
    let d = -self.d;

    DualQuaternion::from_parts(p, d)
  }
}
