use bevy::prelude::*;
use bevy::prelude::shape::UVSphere;
use bevy_aabb_instancing::{
    Cuboid, CuboidMaterial, CuboidMaterialId, CuboidMaterialMap, Cuboids,
    VertexPullingRenderPlugin, COLOR_MODE_SCALAR_HUE,
};
use crate::GRAVITY;

#[derive(Component, Debug, Clone, Copy)]
pub struct Particle {
    pub velocity: Vec3,
    pub density: f32,
    pub pressure_force: Vec3,
    // If particle isn't assigned to a cell, this will be (usize::MAX, usize::MAX, usize::MAX).
    // Order is (layer, row, col)
    pub particle_grid_index: (usize, usize, usize),
    pub predicted_position: Vec3,
    pub position: Vec3,
    pub cuboid: Cuboid,
    // Index of the particle in the vector of particles in the particle grid. Used for checking equality with another particle.
    pub vec_index: usize,
}

#[derive(Component, Debug)]
pub struct ParticleGrid {
    pub cell_width: f32,
    pub num_rows: usize,
    pub num_cols: usize,
    pub num_layers: usize,
    // 3D matrix of lists of particles
    pub particles: Vec<Vec<Vec<Vec<Particle>>>>,
    pub window_width: f32,
    pub window_height: f32,
    pub boundary_dimensions: Vec3,
}

impl Particle {
    pub fn new(cuboid: Cuboid, position: Vec3) -> Self {
        Self {
            velocity: Vec3::new(0.0, 0.0, 0.0),
            density: 0.0,
            pressure_force: Vec3::new(0.0, 0.0, 0.0),
            particle_grid_index: (usize::MAX, usize::MAX, usize::MAX),
            predicted_position: Vec3::new(0.0, 0.0, 0.0),
            position: position,
            cuboid: cuboid,
            vec_index: 0,
        }
    }

    pub fn clone(&self) -> Self {
        Self {
            velocity: self.velocity,
            density: self.density,
            pressure_force: self.pressure_force,
            particle_grid_index: self.particle_grid_index,
            predicted_position: self.predicted_position,
            position: self.position,
            cuboid: self.cuboid.clone(),
            vec_index: self.vec_index,
        }
    }
}

pub fn spawn_particles_grid(
    num_particles: usize,
    particle_radius: f32,
    commands: &mut Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut material_map: ResMut<CuboidMaterialMap>,
    cell_width: f32, num_layers: usize, num_rows: usize, num_cols: usize, window_width: f32, window_height: f32, boundary_dimensions: Vec3
) {
    // Spawn particles in a cube layout.
    // Cube width is limited by the smallest dimension of the boundary - we want particles to fill boundary as much as possible
    // as this looks better.
    // todo: spawn particles to fit the boundary better, rather than always being a cube
    let min_boundary_side_length = f32::min(boundary_dimensions.x, f32::min(boundary_dimensions.y, boundary_dimensions.z));
    let particles_per_side: i32 = (num_particles as f32).powf(1.0 / 3.0).ceil() as i32;
    let spacing: f32 = min_boundary_side_length / particles_per_side as f32;

    let material_id = material_map.push(CuboidMaterial {
        color_mode: COLOR_MODE_SCALAR_HUE,
        ..default()
    });

    let mut instances: Vec<Cuboid> = Vec::with_capacity(num_particles as usize);
    let mut particle_grid = ParticleGrid::new(cell_width, num_layers, num_rows, num_cols, window_width, window_height, boundary_dimensions);
    let mut particles: Vec<Vec<Vec<Vec<Particle>>>> = vec![vec![vec![vec![]; num_cols]; num_rows]; num_layers];
    let mut counter = 0;

    for l in 0..particles_per_side {
        for r in 0..particles_per_side {
            for c in 0..particles_per_side {
                counter += 1;

                let x = ((c - particles_per_side / 2) as f32 + 0.5) as f32 * spacing;
                let y = ((l - particles_per_side / 2) as f32 + 0.5) as f32 * spacing;
                let z = ((r - particles_per_side / 2) as f32 + 0.5) as f32 * spacing;

                let (layer, row, col) = position_to_grid_index(&x, &y, &z, boundary_dimensions, cell_width, num_layers, num_rows, num_cols);
                let scalar_color = u32::from_le_bytes((l * r * c).to_le_bytes());
                let mut cuboid = Cuboid::new(Vec3::new(x - (particle_radius / 2.0), y - (particle_radius / 2.0), z - (particle_radius / 2.0)), Vec3::new(x + (particle_radius / 2.0), y + (particle_radius / 2.0), z + (particle_radius / 2.0)), scalar_color);
                cuboid.set_depth_bias(0);
                instances.push(cuboid);

                let mut particle = Particle::new(cuboid, Vec3::new(x, y, z));
                particle.particle_grid_index = (layer, row, col);
                particle.vec_index = particles[layer as usize][row as usize][col as usize].len();
                particles[layer as usize][row as usize][col as usize].push(particle);
            }

        }
    }
    info!("counter: {}, instances len: {}", counter, instances.len());
    let mut num_particles = 0;

    for l in 0..num_layers {
        for r in 0..num_rows {
            for c in 0..num_cols {
                num_particles += particles[l][r][c].len();
            }
        }
    }

    info!("num_particles: {}", num_particles);

    particle_grid.particles = particles;


    let cuboids = Cuboids::new(instances);

    let aabb = cuboids.aabb();
    commands
        .spawn((SpatialBundle::default(), particle_grid))
        .insert((cuboids, aabb, material_id));
}

impl ParticleGrid {
    pub fn new(cell_width: f32, num_layers: usize, num_rows: usize, num_cols: usize, window_width: f32, window_height: f32, boundary_dimensions: Vec3) -> Self {
        let mut particles = Vec::new();

        for _ in 0..num_layers {
            let mut layer = Vec::new();

            for _ in 0..num_rows {
                let mut row = Vec::new();

                for _ in 0..num_cols {
                        row.push(Vec::new());
                }

                layer.push(row);
            }
            particles.push(layer);
        }

        Self {
            cell_width,
            num_layers: num_layers as usize,
            num_rows: num_rows as usize,
            num_cols: num_cols as usize,
            particles,
            window_width: window_width,
            window_height: window_height,
            boundary_dimensions: boundary_dimensions,
        }
    }

    // Update the grid with the current positions of the particles, and assign each particle to a cell using predicted positions.
    pub fn update_particle_cells_and_predicted_positions(&mut self, boundary_dimensions: Vec3, particle_radius: f32, time: &Res<Time>,) {
        let mut new_particles: Vec<Vec<Vec<Vec<Particle>>>> = vec![vec![vec![vec![]; self.num_cols]; self.num_rows]; self.num_layers];

        for l in 0..self.num_layers {
            for r in 0..self.num_rows {
                for c in 0..self.num_cols {
                    for particle in self.particles[l][r][c].iter_mut() {
                        // assign particle to a cell
                        let (layer, row, col) = position_to_grid_index(&particle.predicted_position.x, &particle.predicted_position.y, &particle.predicted_position.z, self.boundary_dimensions, self.cell_width, self.num_layers, self.num_rows, self.num_cols);
                        particle.particle_grid_index = (layer, row, col);

                        // update predicted position
                        particle.predicted_position.x = f32::max(particle.position.x + particle.velocity.x * time.delta_seconds(), -boundary_dimensions.x / 2.0 + particle_radius);
                        particle.predicted_position.x = f32::min(particle.predicted_position.x, boundary_dimensions.x / 2.0 - particle_radius);
                        particle.predicted_position.y = f32::max(particle.position.y + (particle.velocity.y - GRAVITY) * time.delta_seconds(), -boundary_dimensions.y / 2.0 + particle_radius);
                        particle.predicted_position.y = f32::min(particle.predicted_position.y, boundary_dimensions.y / 2.0 - particle_radius);
                        particle.predicted_position.z = f32::max(particle.position.z + particle.velocity.z * time.delta_seconds(), -boundary_dimensions.z / 2.0 + particle_radius);
                        particle.predicted_position.z = f32::min(particle.predicted_position.z, boundary_dimensions.z / 2.0 - particle_radius);

                        particle.vec_index = new_particles[layer][row][col].len();
                        new_particles[layer][row][col].push(*particle);
                    }
                }
            }
        }

        self.particles = new_particles;
    }

    // Returns the particles in cells within the mouse's sphere of influence.
    // Note: no distance checks are performed.
    pub fn get_particles_adjacent_to_mouse(&self, mouse_pos: Vec3, mouse_radius: &f32) -> Vec<Particle> {
        let (layer, row, col) = position_to_grid_index(&mouse_pos.x, &mouse_pos.y, &0.0, self.boundary_dimensions, self.cell_width, self.num_layers, self.num_rows, self.num_cols);
        let mut neighbours = Vec::new();
        let num_cells_in_radius = (mouse_radius / self.cell_width).ceil() as i32;

        for l in -num_cells_in_radius..=num_cells_in_radius {
            for r in -num_cells_in_radius..=num_cells_in_radius {
                for c in -num_cells_in_radius..=num_cells_in_radius {
                    let new_layer = (layer as i32 + l) as usize;
                    let new_row = (row as i32 + r) as usize;
                    let new_col = (col as i32 + c) as usize;

                    if self.is_row_col_valid(new_layer, new_row, new_col) {
                        neighbours.extend(self.particles[new_layer][new_row][new_col].iter());
                    }
                }
            }
        }

        neighbours
    }

    // Returns the particles in neighbouring cells.
    // Note: The particle itself is not included in the result, and no distance checks are performed.
    pub fn get_particle_neighbours(&self, particle: &Particle, layer: usize, row: usize, col: usize) -> Vec<Particle> {
        let mut neighbours = Vec::new();

        for l in -1..=1 {
            for r in -1..=1 {
                for c in -1..=1 {
                    let new_layer = (layer as i32 + l) as usize;
                    let new_row = (row as i32 + r) as usize;
                    let new_col = (col as i32 + c) as usize;

                    if self.is_row_col_valid(new_layer, new_row, new_col) {
                        // todo: perform distance check here
                        neighbours.extend(self.particles[new_layer][new_row][new_col].iter().filter(|&val| val.vec_index != particle.vec_index));
                    }
                }
            }
        }

        neighbours
    }

    pub fn is_row_col_valid(&self, layer: usize, row: usize, col: usize) -> bool {
        row < self.num_rows && col < self.num_cols && layer < self.num_layers
    }
}

// Returns the row and col index that a particle should be in
fn position_to_grid_index(x: &f32, y: &f32, z: &f32, boundary_dimensions: Vec3, cell_width: f32, num_layers: usize, num_rows: usize, num_cols: usize) -> (usize, usize, usize) {
    let layer = (((y + (boundary_dimensions.y / 2.0)) / cell_width) as usize).clamp(0, num_layers - 1);
    let row = (((z + (boundary_dimensions.z / 2.0)) / cell_width) as usize).clamp(0, num_rows - 1);
    let col = (((x + (boundary_dimensions.x / 2.0)) / cell_width) as usize).clamp(0, num_cols - 1);

    (layer, row, col)
}

pub fn update_scalar_hue_options(time: Res<Time>, mut material_map: ResMut<CuboidMaterialMap>) {
    let material = material_map.get_mut(CuboidMaterialId(1));
    let tv = 1000.0 * (time.elapsed_seconds().sin() + 1.0);
    material.scalar_hue.max_visible = tv;
    material.scalar_hue.clamp_max = tv;
}
