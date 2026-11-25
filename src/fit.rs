use crate::model::Observations;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

// Euclidean projection onto nonnegative weights that sum to one.
fn project(values: &[f64]) -> Vec<f64> {
    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(|a, b| b.total_cmp(a));
    let mut sum = 0.0;
    let mut threshold = 0.0;
    for (i, &v) in sorted.iter().enumerate() {
        sum += v;
        let candidate = (sum - 1.0) / (i + 1) as f64;
        if v <= candidate {
            break;
        }
        threshold = candidate;
    }
    let mut weights: Vec<_> = values.iter().map(|v| (v - threshold).max(0.0)).collect();
    let total: f64 = weights.iter().sum();
    weights.iter_mut().for_each(|x| *x /= total);
    weights
}

impl Observations {
    // Average negative log likelihood per matched marker. On the simplex
    // the other allele has probability 1-p, so one dot product suffices.
    fn evaluate(&self, q: &[f64], mut gradient: Option<&mut [f64]>) -> f64 {
        let mut loss = 0.0;
        for (freq, &[a, b]) in self.frequencies.chunks_exact(q.len()).zip(&self.counts) {
            let p = dot(freq, q);
            let (a, b) = (f64::from(a), f64::from(b));
            loss -= a * p.ln() + b * (-p).ln_1p();
            if let Some(g) = gradient.as_deref_mut() {
                let slope = b / (1.0 - p) - a / p;
                for (g, f) in g.iter_mut().zip(freq) {
                    *g += slope * f;
                }
            }
        }
        let n = self.counts.len() as f64;
        if let Some(g) = gradient {
            g.iter_mut().for_each(|v| *v /= n);
        }
        loss / n
    }
}

pub fn fit(data: &Observations, tolerance: f64) -> Result<Vec<f64>, &'static str> {
    let k = data.components;
    let mut q = vec![1.0 / k as f64; k];
    let mut gradient = vec![0.0; k];
    let mut loss = data.evaluate(&q, Some(&mut gradient));
    let mut step = 1.0;
    for _ in 0..2000 {
        // Convexity bounds the remaining total log-likelihood improvement
        // by n * <gradient, q - e_min>, including optima on the boundary.
        let min = gradient.iter().copied().fold(f64::INFINITY, f64::min);
        let gap = dot(&gradient, &q) - min;
        if gap * data.counts.len() as f64 <= tolerance {
            return Ok(q);
        }
        let proposed: Vec<_> = q
            .iter()
            .zip(&gradient)
            .map(|(v, g)| v - step * (g - min))
            .collect();
        let projected = project(&proposed);
        let direction: Vec<_> = projected
            .iter()
            .zip(&q)
            .map(|(new, old)| new - old)
            .collect();
        let slope = dot(&gradient, &direction);
        let mut scale = 1.0;
        let next = loop {
            let trial: Vec<_> = q
                .iter()
                .zip(&direction)
                .map(|(v, d)| v + scale * d)
                .collect();
            let trial_loss = data.evaluate(&trial, None);
            if trial_loss <= loss + 1e-4 * scale * slope + 1e-13 * loss.abs() {
                break trial;
            }
            scale *= 0.5;
            if scale < 1e-14 {
                return Err("optimizer line search stalled");
            }
        };
        let mut next_gradient = vec![0.0; k];
        loss = data.evaluate(&next, Some(&mut next_gradient));
        let mut squared_distance = 0.0;
        let mut curvature = 0.0;
        for i in 0..k {
            let change = next[i] - q[i];
            squared_distance += change * change;
            curvature += change * (next_gradient[i] - gradient[i]);
        }
        if curvature > 0.0 {
            step = (squared_distance / curvature).clamp(1e-8, 1e5);
        }
        gradient = next_gradient;
        q = next;
    }
    Err("optimizer did not converge within 2000 iterations")
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn matches_analytic_interior_and_boundary_optima() {
        let mut data = Observations {
            components: 2,
            total_markers: 2,
            frequencies: vec![0.9, 0.1, 0.9, 0.1],
            counts: vec![[2, 0], [1, 1]],
        };
        // Three of four calls match the F allele, so 0.1 + 0.8q = 0.75.
        let q = fit(&data, 1e-9).unwrap();
        assert!((q[0] - 0.8125).abs() < 1e-7);
        data.counts[1] = [2, 0];
        assert_eq!(fit(&data, 1e-9).unwrap(), [1.0, 0.0]);
    }
}
