use bevy::{prelude::*, sprite::MaterialMesh2dBundle};

#[derive(Component, Debug, PartialEq)]
pub struct Particle {
    pub velocity: Vec2,
    pub density: f32,
    pub pressure_force: Vec2,
    // if particle isn't assigned to a cell, this will be (usize::MAX, usize::MAX)
    pub particle_grid_index: (usize, usize),
    // pub predicted_position: Vec2,
}

#[derive(Component, Debug)]
pub struct ParticleGrid {
    pub cell_width: f32,
    pub num_rows: usize,
    pub num_cols: usize,
    pub particles: Vec<Vec<Vec<Entity>>>,
    pub window_width: f32,
    pub window_height: f32,
}

impl Particle {
    pub fn new() -> Self {
        Self {
            velocity: Vec2::new(0.0, 0.0),
            density: 0.0,
            pressure_force: Vec2::new(0.0, 0.0),
            particle_grid_index: (usize::MAX, usize::MAX),
            predicted_position: Vec2::new(0.0, 0.0),
        }
    }
}

impl ParticleGrid {
    pub fn new(cell_width: f32, num_rows: usize, num_cols: usize, window: &Window) -> Self {
        let mut particles = Vec::new();

        for _ in 0..num_rows {
            let mut row = Vec::new();

            for _ in 0..num_cols {
                row.push(Vec::new());
            }
            particles.push(row);
        }

        Self {
            cell_width,
            num_rows: num_rows as usize,
            num_cols: num_cols as usize,
            particles,
            window_width: window.width(),
            window_height: window.height(),
        }
    }

    pub fn spawn_particles(
        &mut self,
        num_particles: usize,
        particle_radius: f32,
        commands: &mut Commands,
        mut meshes: ResMut<Assets<Mesh>>,
        mut materials: ResMut<Assets<ColorMaterial>>,
    ) {
        // spawn particles in a grid
        let particles_per_row: i32 = (num_particles as f32).sqrt() as i32;
        let particles_per_col: i32 = (num_particles as i32 - 1) / particles_per_row + 1;
        let spacing: f32 = particle_radius * 6.0;

        for i in 0..num_particles {
            let x = ((i as i32 % particles_per_row - particles_per_row / 2) as f32 + 0.5) as f32 * spacing;
            let y = ((i as i32 / particles_per_row - particles_per_col / 2) as f32 + 0.5) as f32 * spacing;
            let (row, col) = self.position_to_grid_index(&x, &y);

            let mut particle = Particle::new();
            particle.particle_grid_index = (row, col);

            let transform = Transform::from_translation(Vec3::new(x, y, 0.0));

            let entity = commands
                .spawn(particle)
                .insert(MaterialMesh2dBundle {
                    mesh: meshes.add(shape::Circle {
                        radius: particle_radius,
                        ..Default::default()
                    }.into()).into(),
                    material: materials.add(ColorMaterial::from(Color::WHITE)),
                    transform: transform,
                    ..Default::default()
                }).id();

            self.particles[row][col].push(entity);
        }
    }

    pub fn update_particle_cells(&mut self, particles: &mut Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>) {
        let particles_clone = self.particles.clone();

        for r in 0..self.num_rows {
            for c in 0..self.num_cols {
                self.particles[r][c].clear();
            }
        }

        for r in 0..self.num_rows {
            for c in 0..self.num_cols {
                for entity in particles_clone[r][c].iter() {
                    let (mut particle, transform, _) = particles.get_mut(*entity).unwrap();

                    // assign particle to a cell
                    let (row, col) = self.position_to_grid_index(&transform.translation.x, &transform.translation.y);
                    particle.particle_grid_index = (row, col);
                    self.particles[row][col].push(*entity);
                }
            }
        }
    }

    // Returns the row and col index that a particle should be in
    pub fn position_to_grid_index(&self, x: &f32, y: &f32) -> (usize, usize) {
        let row = ((y + (self.window_height / 2.0)) / self.cell_width) as usize;
        let col = ((x + (self.window_width / 2.0)) / self.cell_width) as usize;

        (row, col)
    }

    pub fn is_row_col_valid(&self, row: usize, col: usize) -> bool {
        row < self.num_rows && col < self.num_cols
    }
}
