// Custom SPH function for viscosity force.
pub fn viscosity_smoothing_kernel(radius: &f32, dst: &f32) -> f32 {
    let value = f32::max(0., radius * radius - dst * dst);
    value * value * value / (radius * radius * radius * radius * radius * radius)
}

// Calculate the relative 'influence' of a particle at a given distance from a point.
pub fn smoothing_kernel(radius: &f32, dst: &f32) -> f32 {
    if dst >= radius {
        return 0.0;
    }

    let volume = (std::f32::consts::PI * radius.powi(4)) / 6.0;
    return (radius - dst) * (radius - dst) / volume;
}

// Calculates gradient of the smoothing kernel at a given distance.
pub fn smoothing_kernel_derivative(radius: &f32, dst: &f32) -> f32 {
    let scale = 12. / (std::f32::consts::PI * radius.powi(4));
    return scale * (dst - radius);
}