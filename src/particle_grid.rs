use bevy::{prelude::*, sprite::MaterialMesh2dBundle};

#[derive(Component, Debug, PartialEq)]
pub struct Particle {
    pub velocity: Vec2,
    pub density: f32,
    pub pressure_force: Vec2,
    // if particle isn't assigned to a cell, this will be (usize::MAX, usize::MAX)
    pub particle_grid_index: (usize, usize),
    pub predicted_position: Vec2,
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

    // Update the grid with the current positions of the particles, and assign each particle to a cell using predicted positions.
    pub fn update_particle_cells(&mut self, particles: &mut Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>) {
        // todo: find a way to do this without cloning
        let particles_clone = self.particles.clone();

        for r in 0..self.num_rows {
            for c in 0..self.num_cols {
                self.particles[r][c].clear();
            }
        }

        for r in 0..self.num_rows {
            for c in 0..self.num_cols {
                for entity in particles_clone[r][c].iter() {
                    let (mut particle, _, _) = particles.get_mut(*entity).unwrap();

                    // assign particle to a cell
                    // let (row, col) = self.position_to_grid_index(&transform.translation.x, &transform.translation.y);
                    let (row, col) = self.position_to_grid_index(&particle.predicted_position.x, &particle.predicted_position.y);
                    particle.particle_grid_index = (row, col);
                    self.particles[row][col].push(*entity);
                }
            }
        }
    }

    // Returns the entities of particles in cells within the mouse's sphere of influence.
    // Note: no distance checks are performed.
    pub fn get_particles_adjacent_to_mouse(&self, mouse_pos: Vec2, mouse_radius: &f32) -> Vec<Entity> {
        let (row, col) = self.position_to_grid_index(&mouse_pos.x, &mouse_pos.y);
        let mut neighbours = Vec::new();
        let num_cells_in_radius = (mouse_radius / self.cell_width).ceil() as i32;

        for r in -num_cells_in_radius..=num_cells_in_radius {
            for c in -num_cells_in_radius..=num_cells_in_radius {
                let new_row = (row as i32 + r) as usize;
                let new_col = (col as i32 + c) as usize;

                if self.is_row_col_valid(new_row, new_col) {
                    neighbours.extend(self.particles[new_row][new_col].iter());
                }
            }
        }

        neighbours
    }

    // Returns the entities of particles in neighbouring cells.
    // Note: The particle itself is not included in the result, and no distance checks are performed.
    pub fn get_particle_neighbours(&self, entity: &Entity, row: usize, col: usize) -> Vec<Entity> {
        let mut neighbours = Vec::new();

        for r in -1..=1 {
            for c in -1..=1 {
                let new_row = (row as i32 + r) as usize;
                let new_col = (col as i32 + c) as usize;

                if self.is_row_col_valid(new_row, new_col) {
                    neighbours.extend(self.particles[new_row][new_col].iter().filter(|&&e| e != *entity));
                }
            }
        }

        neighbours
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
