use std::ops::Deref;

use crate::{Interpolation, Interpolator};

pub struct DynInterpolation<T> {
    #[allow(clippy::type_complexity)]
    build: Box<dyn (Fn(&[T], &[T]) -> Result<Box<dyn Interpolator<T>>, crate::InterpolationError>)>,
    name: Box<str>,
    min_size: usize,
}

impl<T> DynInterpolation<T> {
    pub fn new<I>(interp: I) -> Self
    where
        I: Interpolation<T> + 'static,
        I::Interpolator: 'static,
    {
        Self {
            name: interp.name().into(),
            min_size: interp.min_size(),
            build: Box::new(move |xa, ya| match interp.build(xa, ya) {
                Ok(interp) => Ok(Box::new(interp)),
                Err(err) => Err(err),
            }),
        }
    }
}

impl<T> Interpolation<T> for DynInterpolation<T> {
    type Interpolator = Box<dyn Interpolator<T>>;

    fn build(&self, xa: &[T], ya: &[T]) -> Result<Self::Interpolator, crate::InterpolationError> {
        (self.build)(xa, ya)
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn min_size(&self) -> usize {
        self.min_size
    }
}

impl<T> Interpolator<T> for Box<dyn Interpolator<T>> {
    fn eval(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval(xa, ya, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv(xa, ya, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_deriv2(xa, ya, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut crate::Accelerator,
    ) -> Result<T, crate::DomainError> {
        self.deref().eval_integ(xa, ya, a, b, acc)
    }
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_dyn_interpolation() {
        let xa = [1.0, 2.0, 3.0];
        let ya = [1.0, 2.0, 3.0];

        let interp: Box<dyn Interpolator<_>> =
            DynInterpolation::new(Cubic).build(&xa, &ya).unwrap();

        fn eval<I: Interpolator<T>, T>(interp: I, xa: &[T], ya: &[T], x: T) {
            let mut accel = Accelerator::new();
            let _ = interp.eval(xa, ya, x, &mut accel).unwrap();
        }

        eval(interp, &xa, &ya, 1.5);
    }
}
