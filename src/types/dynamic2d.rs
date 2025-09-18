use std::ops::Deref;

use crate::Interp2dType;
use crate::Interpolation2d;
use crate::InterpolationError;

/// Representation of a 2d Interpolation Type that is not known in compile-time.
pub struct DynInterp2dType<T> {
    #[allow(clippy::type_complexity)]
    build:
        Box<dyn (Fn(&[T], &[T], &[T]) -> Result<Box<dyn Interpolation2d<T>>, InterpolationError>)>,
    name: Box<str>,
    min_size: usize,
}

impl<T> DynInterp2dType<T> {
    pub fn new<I>(interp: I) -> Self
    where
        I: Interp2dType<T> + 'static,
        I::Interpolation2d: 'static,
    {
        Self {
            name: interp.name().into(),
            min_size: interp.min_size(),
            build: Box::new(move |xa, ya, za| match interp.build(xa, ya, za) {
                Ok(interp) => Ok(Box::new(interp)),
                Err(err) => Err(err),
            }),
        }
    }
}

impl<T> Interp2dType<T> for DynInterp2dType<T> {
    type Interpolation2d = Box<dyn Interpolation2d<T>>;

    fn build(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
    ) -> Result<Self::Interpolation2d, InterpolationError> {
        (self.build)(xa, ya, za)
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn min_size(&self) -> usize {
        self.min_size
    }
}

impl<T> Interpolation2d<T> for Box<dyn Interpolation2d<T>> {
    fn eval_extrap(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut crate::Accelerator,
        yacc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_extrap(xa, ya, za, x, y, xacc, yacc)
    }

    fn eval_deriv_x(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut crate::Accelerator,
        yacc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv_x(xa, ya, za, x, y, xacc, yacc)
    }

    fn eval_deriv_y(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut crate::Accelerator,
        yacc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv_y(xa, ya, za, x, y, xacc, yacc)
    }

    fn eval_deriv_xx(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut crate::Accelerator,
        yacc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv_xx(xa, ya, za, x, y, xacc, yacc)
    }

    fn eval_deriv_yy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut crate::Accelerator,
        yacc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv_yy(xa, ya, za, x, y, xacc, yacc)
    }

    fn eval_deriv_xy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut crate::Accelerator,
        yacc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv_xy(xa, ya, za, x, y, xacc, yacc)
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn test_dyn_interp2d_type() {
        let xa = [0.0, 1.0, 2.0, 3.0];
        let ya = [0.0, 1.0, 2.0, 3.0];
        #[rustfmt::skip]
        let za = [
            1.0, 1.1, 1.2, 1.3,
            1.1, 1.2, 1.3, 1.4,
            1.2, 1.3, 1.4, 1.5,
            1.3, 1.4, 1.5, 1.6,
        ];

        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let interp: Box<dyn Interpolation2d<_>> =
            DynInterp2dType::new(Bicubic).build(&xa, &ya, &za).unwrap();

        let _ = interp
            .eval(&xa, &ya, &za, 1.0, 1.0, &mut xacc, &mut yacc)
            .unwrap();
    }
}
