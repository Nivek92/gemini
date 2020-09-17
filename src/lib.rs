#![feature(clamp)]
mod geometry {

    use nalgebra::base::dimension::{U3, U4, U6, U8};
    use nalgebra::base::{VectorN, MatrixN};
    use nalgebra::geometry::Quaternion;
    use std::ops::{Add, Mul, Sub, Div};

    static THRESHOLD: f64 = 1e-12;

    pub type Vector3 = VectorN<f64, U3>;
    pub type Vector4 = VectorN<f64, U4>;
    pub type Vector6 = VectorN<f64, U6>;
    pub type Vector8 = VectorN<f64, U8>;
    pub type Matrix4 = MatrixN<f64, U4>;
    pub type Matrix8 = MatrixN<f64, U8>;

    pub struct DualQuaternion {
        pub p: Quaternion<f64>,
        pub d: Quaternion<f64>,
    }

    impl DualQuaternion {

        pub fn new(q0: f64, q1: f64, q2: f64, q3: f64, q4: f64, q5: f64, q6: f64, q7: f64) -> Self {
            DualQuaternion {
                p: Quaternion::new(q0, q1, q2, q3),
                d: Quaternion::new(q4, q5, q6, q7),
            }
        }

        pub fn from(v: Vector8) -> Self {
            DualQuaternion {
                p: Quaternion::new(v[0], v[1], v[2], v[3]),
                d: Quaternion::new(v[4], v[5], v[6], v[7]),
            }
        }

        pub fn p(&self) -> DualQuaternion {
            DualQuaternion::new(self.p[3], self.p[0], self.p[1], self.p[2], 0.0, 0.0, 0.0, 0.0)
        }

        pub fn primary(&self) -> DualQuaternion {
            self.p()
        }

        pub fn d(&self) -> DualQuaternion {
            DualQuaternion::new(self.d[3], self.d[0], self.d[1], self.d[2], 0.0, 0.0, 0.0, 0.0)
        }

        pub fn dual(&self) -> DualQuaternion {
            self.d()
        }

        pub fn re(&self) -> DualQuaternion {
            DualQuaternion::new(
                self.p[3],
                0.0,
                0.0,
                0.0,
                self.d[3],
                0.0,
                0.0,
                0.0,
            )
        }

        pub fn real(&self) -> DualQuaternion {
            self.re()
        }

        pub fn im(&self) -> DualQuaternion {
            DualQuaternion::new(
                0.0,
                self.p[0],
                self.p[1],
                self.p[2],
                0.0,
                self.d[0],
                self.d[1],
                self.d[2],
            )
        }

        pub fn imaginary(&self) -> DualQuaternion {
            self.im()
        }

        pub fn conj(&self) -> DualQuaternion {
            DualQuaternion::new(
                self.p[3], -self.p[0], -self.p[1], -self.p[2], self.d[3], -self.d[0], -self.d[1],
                -self.d[2],
            )
        }

        pub fn conjugate(&self) -> DualQuaternion {
            self.conj()
        }

        pub fn norm(self) -> DualQuaternion {
            if self.is_pure() {
                return DualQuaternion::from(Vector8::zeros());
            }

            let mut values = Vector8::zeros();
            let norm = self.conj() * self;
            let _p = norm.p;
            let _d = norm.d;


            values[0] = f64::sqrt(_p[3]);
            values[4] = _d[3] / 2.0 * _p[3];

            for i in 0..4 {
                if _p[i] < THRESHOLD {
                    values[i] = 0.0;
                }
            }

            for i in 0..4 {
                if _d[i] < THRESHOLD {
                    values[4 + i] = 0.0;
                }
            }

            DualQuaternion::from(Vector8::zeros())
        }

        pub fn inv(self) -> DualQuaternion{
            let norm = self * self.conj();

            self.conj() * DualQuaternion::new( 1.0 / norm.p[3], 0.0, 0.0, 0.0, -norm.d[3] / (norm.p[3] * norm.p[3]), 0.0, 0.0, 0.0)
        }

        pub fn translation(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The translation operation is defined only for unit dual quaterions.")
            }

            2. * (self.d() * self.p().conjugate())
        }

        pub fn rotation(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The rotation operation is defined only for unit dual quaterions.")
            }

            self.p()
        }

        pub fn rotation_axis(&self) -> DualQuaternion {
            // if !self.is_unit() {
            //     panic!("The rotation_axis operation is defined only for unit dual quaterions.")
            // }

            // let phi = rotation_angle() / 2.0;
            // if(phi == 0.0) {
            //     return DualQuaternion::new(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0); // DQ(0,0,0,1)
            // }  
            // else
            // {
            //     //rotation axis calculation
            //     let rot_axis = self.p();
            //     rot_axis = ( rot_axis.im() * (1/sin(phi)) );

            //     return rot_axis;
            // }
            unimplemented!
        }

        pub fn rotation_angle(&self) -> f64 {
            if !self.is_unit() {
                panic!("The rotation_angle operation is defined only for unit dual quaterions.")
            }

            2. * f64::acos(self.p[3].clamp(-1.0, 1.0))
        }

        pub fn log(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The log operation is defined only for unit dual quaterions.")
            }

            let q1: DualQuaternion = (0.5 * self.rotation_angle()) * self.rotation_axis();
            let q2: DualQuaternion = 0.5 * self.translation();

            let p: Quaternion<f64> = q1.p;
            let d: Quaternion<f64> = q2.p;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }

        pub fn exp(&self) -> DualQuaternion {
            if !self.is_pure() {
                panic!("The exponential operation is defined only for pure dual quaterions.")
            }

            let mut prim = self.p();
             let phi = prim.p.norm();

            if phi != 0.0 {
                prim = f64::cos(phi) + (f64::sin(phi) / phi) * prim;
            } else {
                prim = DualQuaternion::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            }

            prim + DualQuaternion::new(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0) * self.d() * prim
        }

        pub fn pow(&self, a: f64) -> DualQuaternion {
            (a * self.log()).exp()
        }

        pub fn tplus(self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The tplus operation is defined only for unit dual quaterions.")
            }

            self * self.p().conj()
        }

        pub fn pinv(&self) -> DualQuaternion {
            if !self.is_unit() {
                panic!("The pinv operation is defined only for unit dual quaterions.")
            }

            (self.conj().tplus() * self.tplus()).conj() * self.conj()
        }

        pub fn ham_add4(&self) -> Matrix4 {
            Matrix4::new(
                self.p[3], -self.p[0], -self.p[1], -self.p[2],
                self.p[0],  self.p[3], -self.p[2],  self.p[1],
                self.p[1],  self.p[2],  self.p[3], -self.p[0],
                self.p[2], -self.p[1],  self.p[0],  self.p[3],
            )
        }

        pub fn ham_subtract4(&self) -> Matrix4 {
            Matrix4::new(
                self.p[3], -self.p[0], -self.p[1], -self.p[2],
                self.p[0],  self.p[3],  self.p[2], -self.p[1],
                self.p[1], -self.p[2],  self.p[3],  self.p[0],
                self.p[2],  self.p[1], -self.p[0],  self.p[3],
            )
        }

        pub fn ham_add8(&self) -> Matrix8 {
            Matrix8::from_row_slice(&[
                self.p[3], -self.p[0], -self.p[1], -self.p[2], 0.0, 0.0, 0.0, 0.0,
                self.p[0], self.p[3], -self.p[2], self.p[1], 0.0, 0.0, 0.0, 0.0,
                self.p[1], self.p[2], self.p[3], -self.p[0], 0.0, 0.0, 0.0, 0.0,
                self.p[2], -self.p[1], self.p[0], self.p[3], 0.0, 0.0, 0.0, 0.0,
                self.d[3], -self.d[0], -self.d[1], -self.d[2], self.p[3], -self.p[0], -self.p[1], -self.p[2],
                self.d[0], self.d[3], -self.d[2], self.d[1], self.p[0], self.p[3], -self.p[2], self.p[1],
                self.d[1], self.d[2], self.d[3], -self.d[0], self.p[1], self.p[2], self.p[3], -self.p[0],
                self.d[2], -self.d[1], self.d[0], self.d[3], self.p[2], -self.p[1], self.p[0], self.p[3],
            ])
        }

        pub fn ham_subtract8(&self) -> Matrix8 {
            Matrix8::from_row_slice(&[
                self.p[3], -self.p[0], -self.p[1], -self.p[2], 0.0, 0.0, 0.0, 0.0,
                self.p[0], self.p[3], self.p[2], -self.p[1], 0.0, 0.0, 0.0, 0.0,
                self.p[1], -self.p[2], self.p[3], self.p[0], 0.0, 0.0, 0.0, 0.0,
                self.p[2], self.p[1], -self.p[0], self.p[3], 0.0, 0.0, 0.0, 0.0,
                self.d[3], -self.d[0], -self.d[1], -self.d[2], self.p[3], -self.p[0], -self.p[1], -self.p[2],
                self.d[0], self.d[3], self.d[2], -self.d[1], self.p[0], self.p[3], self.p[2], -self.p[1],
                self.d[1], -self.d[2], self.d[3], self.d[0], self.p[1], -self.p[2], self.p[3], self.p[0],
                self.d[2], self.d[1], -self.d[0], self.d[3], self.p[2], self.p[1], -self.p[0], self.p[3],
            ])
        }

        pub fn vec3(&self) -> Vector3 {
            Vector3::new(self.p[0], self.p[1], self.p[2])
        }

        pub fn vec4(&self) -> Vector4 {
            Vector4::new(self.p[3], self.p[0], self.p[1], self.p[2])
        }

        pub fn vec6(&self) -> Vector6 {
            Vector6::new(self.p[0], self.p[1], self.p[2], self.d[0], self.d[1], self.d[2])
        }

        pub fn vec8(&self) -> Vector8 {
            Vector8::from_row_slice(&[self.p[3], self.p[0], self.p[1], self.p[2], self.d[3], self.d[0], self.d[1], self.d[2]])
        }

        pub fn generalized_jacobian(&self) -> Matrix8 {
            Matrix8::from_row_slice(&[
                self.d[3],  self.d[0],  self.d[1],  self.d[2],  self.p[3],  self.p[0],  self.p[1],  self.p[2],
                self.d[0], -self.d[3],  self.d[2], -self.d[1], -self.p[0],  self.p[3], -self.p[2],  self.p[1],
                self.d[1], -self.d[2], -self.d[3],  self.d[0], -self.p[1],  self.p[2],  self.p[3], -self.p[0],
                self.d[2],  self.d[1], -self.d[0], -self.d[3], -self.p[2], -self.p[1],  self.p[0],  self.p[3],
                self.p[3],  self.p[0],  self.p[1],  self.p[2], 0.0, 0.0, 0.0, 0.0,
                -self.p[0],  self.p[3], -self.p[2],  self.p[1], 0.0, 0.0, 0.0, 0.0,
                -self.p[1],  self.p[2],  self.p[3], -self.p[0], 0.0, 0.0, 0.0, 0.0,
                -self.p[2], -self.p[1],  self.p[0],  self.p[3], 0.0, 0.0, 0.0, 0.0,
            ])
        }

        pub fn normalize(self) -> DualQuaternion {
            self * self.norm().inv()
        }

        pub fn sharp(&self) -> DualQuaternion {
            self.p() - DualQuaternion::new(0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0) * self.d()
        }

        pub fn adj() {}

        pub fn adj_sharp() {}

        pub fn to_unit(&self) {
            
        }

        pub fn is_unit(&self) -> bool {
            unimplemented!
        }

        pub fn is_pure(&self) -> bool {
            let r = self.re();
            r.p[3] == 0.0 && r.d[3] == 0.0
        }
        pub fn is_real(&self) -> bool {
            let r = self.re();
            // TODO: replace with methods on individual quaternions _p.is_real() && _d.is_real()
            r.p[0] == 0.0
                && r.p[1] == 0.0
                && r.p[2] == 0.0
                && r.d[0] == 0.0
                && r.d[1] == 0.0
                && r.d[2] == 0.0
        }

        pub fn is_real_number(&self) -> bool {
            self.is_real() && self.is_quaternion()
        }

        pub fn is_quaternion(&self) -> bool {
            self.d == Quaternion::default()
        }

        pub fn is_pure_quaternion(&self) -> bool {
            self.is_pure() && self.is_quaternion()
        }

        pub fn is_line(&self) -> bool {
            self.is_unit() && self.is_pure()
        }

        pub fn is_plane(&self) -> bool {
            let is_real = self.d[0] == 0.0 && self.d[1] == 0.0 && self.d[2] == 0.0; // TODO: add method to nalgebra library and call it instead
            self.is_unit() && is_real
        }
    }

    // Addition: DualQuaternion
    impl Add<DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn add(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p + q2.p;
            let d = q1.d + q2.d;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }
    }

    // Addition: DualQuaternion Reference
    impl Add<&DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn add(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p + q2.p;
            let d = q1.d + q2.d;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }
    }

    // Addition: RHS scalar
    impl Add<f64> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn add(self, rhs: f64) -> DualQuaternion {
            self + DualQuaternion::new(rhs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        }
    }

    // Addition: LHS scalar
    impl Add<DualQuaternion> for f64 {
        type Output = DualQuaternion;

        #[inline]
        fn add(self, rhs: DualQuaternion) -> DualQuaternion {
            DualQuaternion::new(self, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) + rhs
        }
    }

    // Subtraction: DualQuaternion
    impl Sub<DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn sub(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p - q2.p;
            let d = q1.d - q2.d;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }
    }

    // Substraction: RHS scalar
    impl Sub<f64> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn sub(self, rhs: f64) -> DualQuaternion {
            self - DualQuaternion::new(rhs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        }
    }

    // Substraction: LHS scalar
    impl Sub<DualQuaternion> for f64 {
        type Output = DualQuaternion;

        #[inline]
        fn sub(self, rhs: DualQuaternion) -> DualQuaternion {
            DualQuaternion::new(self, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) - rhs
        }
    }

    // Subtraction: DualQuaternion Reference
    impl Sub<&DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        fn sub(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p - q2.p;
            let d = q1.d - q2.d;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }
    }

    // Multiplication: DualQuaternion
    impl Mul<DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, q2: DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p * q1.p;
            let d = q1.p * q2.d + q1.d * q2.p;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }
    }

    // Multiplaction: DualQuaternion Reference
    impl Mul<&DualQuaternion> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, q2: &DualQuaternion) -> DualQuaternion {
            let q1 = &self;
            let p = q1.p * q1.p;
            let d = q1.p * q2.d + q1.d * q2.p;

            DualQuaternion::new(p[3], p[0], p[1], p[2], d[3], d[0], d[1], d[2])
        }
    }

    // Multiplaction: RHS scalar
    impl Mul<f64> for DualQuaternion {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, rhs: f64) -> DualQuaternion {
            self * DualQuaternion::new(rhs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        }
    }

    // Multiplaction: LHS scalar
    impl Mul<DualQuaternion> for f64 {
        type Output = DualQuaternion;

        #[inline]
        fn mul(self, rhs: DualQuaternion) -> DualQuaternion {
            DualQuaternion::new(self, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * rhs
        }
    }
}
