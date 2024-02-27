mod particle_grid;
// used for better error messages in web assembly
extern crate console_error_panic_hook;

use bevy::prelude::*;
use bevy::ecs::query::BatchingStrategy;
use bevy::window::{PrimaryWindow, PresentMode};
use bevy::utils::Instant;
use rand::Rng;
use particle_grid::{Particle, ParticleGrid};
use rayon::prelude::*;

/// We will store the world position of the mouse cursor here.
#[derive(Resource, Default)]
struct MyWorldCoords(Vec2);

/// Used to help identify our main camera
#[derive(Component)]
struct MainCamera;

#[derive(Component)]
struct Timer {
    start_time: Instant,
    // how many iterations to sum time for before printing the average
    num_iters: u32,
}

const NUM_PARTICLES: usize = 1052;
const VISUALIZE_COLOR_BASED_ON: &str = "velocity"; // density or velocity
const PARTICLE_RADIUS: f32 = 3.;
const RESPOND_TO_MOUSE: bool = true;
const MOUSE_RADIUS: f32 = 150.0;
// how much force the mouse applies to the particles
const MOUSE_PRESSURE_MULTIPLIER: f32 = 10.0;
const GRAVITY: f32 = 0.;
const WINDOW_WIDTH: f32 = 800.0;
const WINDOW_HEIGHT: f32 = 600.0;
const COLLISION_DAMPING: f32 = 0.6;
const PARTICLE_MASS: f32 = 1.0;
const SMOOTHING_RADIUS: f32 = 50.;
const TARGET_DENSITY: f32 = 0.0018;
const PRESSURE_MULTIPLIER: f32 = 2500.;
const VISCOSITY_STRENGTH: f32 = 0.007;
const RUN_STOPWATCH: bool = false;
// how many times to run the simulation before printing the average time to update positions
const NUM_ITERATIONS_BEFORE_PRINT_STOPWATCH: u32 = 300;
// How much the window edge repels particles.
// X value is for left/right walls, Y value is for top/bottom walls.
const BOUNDARY_PRESSURE_MULTIPLIER: Vec2 = Vec2::new(0.01, 0.01);
const BUILD_FOR_WEB: bool = true;
// const BOUNDARY_PRESSURE_MULTIPLIER: Vec2 = Vec2::new(0.02, 0.5);

fn main() {
    console_error_panic_hook::set_once();
    
    let mut app = App::new();

    app
        .init_resource::<MyWorldCoords>()
        .insert_resource(ClearColor(Color::rgb(0.0, 0.0, 0.0)))
        .add_systems(Startup, setup)
        .add_systems(FixedUpdate, (update_particle_positions, update_densities));

    if BUILD_FOR_WEB == false {
        app.add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Fluid Simulation".into(),
                resolution: (WINDOW_WIDTH, WINDOW_HEIGHT).into(),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }));
    } else {
        // Tell the app the ID of the canvas element we want to use.
        // This ID needs to match the ID of the canvas element in the HTML file.
        app.add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                // provide the ID selector string here
                canvas: Some("#fluid-simulation-canvas".into()),
                title: "Fluid Simulation".into(),
                resolution: (WINDOW_WIDTH, WINDOW_HEIGHT).into(),
                ..default()
            }),
            ..default()
        }));
    }

    app.run();
}

fn setup(
    mut commands: Commands,
    meshes: ResMut<Assets<Mesh>>,
    materials: ResMut<Assets<ColorMaterial>>,
) {
    commands.spawn((Camera2dBundle::default(), MainCamera, Timer { start_time: Instant::now(), num_iters: NUM_ITERATIONS_BEFORE_PRINT_STOPWATCH }));
    spawn_particles(commands, meshes, materials);
}

fn spawn_particles(
    mut commands: Commands,
    meshes: ResMut<Assets<Mesh>>,
    materials: ResMut<Assets<ColorMaterial>>,
) {
    // create spatial partitioning grid
    let mut particle_grid = ParticleGrid::new(SMOOTHING_RADIUS, (WINDOW_HEIGHT / SMOOTHING_RADIUS).ceil() as usize, (WINDOW_WIDTH / SMOOTHING_RADIUS).ceil() as usize, WINDOW_WIDTH, WINDOW_HEIGHT);
    particle_grid.spawn_particles(NUM_PARTICLES, PARTICLE_RADIUS, &mut commands, meshes, materials);
    commands.spawn(particle_grid);
}

// Updates the position of the particles based on their density and pressure.
fn update_particle_positions(
    mut particle_grid_query: Query<&mut ParticleGrid>,
    mut particles: Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>,
    mut mycoords: ResMut<MyWorldCoords>, 
    mut timer_query: Query<&mut Timer>,
    q_window: Query<&Window, With<PrimaryWindow>>, 
    q_camera: Query<(&Camera, &GlobalTransform), With<MainCamera>>, 
    time: Res<Time>,
) {
    if RUN_STOPWATCH {
        let mut timer = timer_query.get_single_mut().unwrap();

        if timer.num_iters == 0 {
            info!("Average time per iteration: {:?}", timer.start_time.elapsed() / NUM_ITERATIONS_BEFORE_PRINT_STOPWATCH);
            timer.start_time = Instant::now();
            timer.num_iters = NUM_ITERATIONS_BEFORE_PRINT_STOPWATCH;
        } else {
            timer.num_iters -= 1;
        }
    }

    // let start_time = Instant::now();
    let window_state = q_window.get_single().unwrap();
    let mut grid = particle_grid_query.single_mut();

    // Calcualte each particle's predicted position to use for improved pressure force calculation.
    // Choose batch size of 32 to limit overhead of ParallelIterator, since getting predicted pos is inexpensive.
    particles
    .par_iter_mut()
    .batching_strategy(BatchingStrategy::fixed(32))
    .for_each(|(mut particle, transform, _)| {
        particle.predicted_position.x = f32::max(transform.translation.x + particle.velocity.x * time.delta_seconds(), -WINDOW_WIDTH / 2.0 + PARTICLE_RADIUS);
        particle.predicted_position.x = f32::min(particle.predicted_position.x, WINDOW_WIDTH / 2.0 - PARTICLE_RADIUS);
        particle.predicted_position.y = f32::max(transform.translation.y + (particle.velocity.y - GRAVITY) * time.delta_seconds(), -WINDOW_HEIGHT / 2.0 + PARTICLE_RADIUS);
        particle.predicted_position.y = f32::min(particle.predicted_position.y, WINDOW_HEIGHT / 2.0 - PARTICLE_RADIUS);
    });

    grid.update_particle_cells(&mut particles);

    if RESPOND_TO_MOUSE {
        let mouse_pos = get_mouse_world_position(&mut mycoords, &window_state, &q_camera);

        if mouse_pos != None {
            // find particles within mouse's radius and apply a force to them
            let neighbours = grid.get_particles_adjacent_to_mouse(mycoords.0, &MOUSE_RADIUS);

            unsafe {
                neighbours.par_iter().for_each(|entity| {
                        // for entity in neighbours.iter() {
                    // This line is the reason we need the "unsafe" above.
                    let (mut particle, transform, _) = particles.get_unchecked(*entity).unwrap();
                    let dir = mycoords.0 - transform.translation.xy();
                    let euclidean_distance = dir.length();
                    
                    if euclidean_distance <= MOUSE_RADIUS {
                        particle.velocity += smoothing_kernel(&MOUSE_RADIUS, &euclidean_distance) * dir * MOUSE_PRESSURE_MULTIPLIER * PRESSURE_MULTIPLIER * time.delta_seconds();
                    }
                });
            }
        }
    }

    update_pressure_forces(&mut particles, &grid);

    particles
    .par_iter_mut()
    .batching_strategy(BatchingStrategy::fixed(32))
    .for_each(|(mut particle, mut transform, _)| {
        // apply gravity
        particle.velocity.y -= GRAVITY * time.delta_seconds();

        if particle.density != 0.0 {
            // F = m * a, so a = F / m
            let pressure_force = particle.pressure_force / particle.density;
            particle.velocity += pressure_force * time.delta_seconds();
        }
        
        transform.translation.x += particle.velocity.x;
        transform.translation.y += particle.velocity.y;

        resolve_collisions(&mut particle, &mut transform);
    });
}

// Returns the mouse position in world coordinates
fn get_mouse_world_position(
    mycoords: &mut ResMut<MyWorldCoords>,
    // query to get the window (so we can read the current cursor position)
    window_state: &Window,
    // query to get camera transform
    q_camera: &Query<(&Camera, &GlobalTransform), With<MainCamera>>,
) -> Option<Vec2> {
    // get the camera info and transform
    // assuming there is exactly one main camera entity, so Query::single() is OK
    let (camera, camera_transform) = q_camera.single();

    // check if the cursor is inside the window and get its position
    // then, ask bevy to convert into world coordinates, and truncate to discard Z
    if let Some(world_position) = window_state.cursor_position()
        .and_then(|cursor| camera.viewport_to_world(camera_transform, cursor))
        .map(|ray| ray.origin.truncate())
    {
        mycoords.0 = world_position;
        return Some(world_position);
    }

    // if the cursor is not inside the window, we don't update the world position
    return None;
}

// Returns the force acting on a particle due to its proximity to the edge of the window.
// For simplicity, assume the radius of influence of the window edge is the same as the smoothing radius.
fn calculate_boundary_force(transform: &Transform) -> Vec2 {
    let mut force = Vec2::new(0.0, 0.0);
    let dist_from_x_boundary = (WINDOW_WIDTH / 2.0) - transform.translation.x.abs();
    let dist_from_y_boundary = (WINDOW_HEIGHT / 2.0) - transform.translation.y.abs();

    if dist_from_x_boundary < SMOOTHING_RADIUS {
        force.x = BOUNDARY_PRESSURE_MULTIPLIER.x * viscosity_smoothing_kernel(&SMOOTHING_RADIUS, &dist_from_x_boundary) * -sign(transform.translation.x);
    }

    if dist_from_y_boundary < SMOOTHING_RADIUS {
        force.y = BOUNDARY_PRESSURE_MULTIPLIER.y * viscosity_smoothing_kernel(&SMOOTHING_RADIUS, &dist_from_y_boundary) * -sign(transform.translation.y);
    }

    force
}

// Calculate the force between all particles to simulate pressure.
fn update_pressure_forces(particles: &mut Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>, particle_grid: &ParticleGrid) {
    // for row in 0..particle_grid.num_rows {
        particle_grid.particles.par_iter().enumerate().for_each(|(row, _)| {

        for col in 0..particle_grid.num_cols {
            for entity in particle_grid.particles[row][col].iter() {
            // particle_grid.particles[row][col].par_iter().for_each(|entity| {
                unsafe {
                    let (mut particle, transform, _) = particles.get_unchecked(*entity).unwrap();
                    let neighbours = particle_grid.get_particle_neighbours(&entity, row, col);
                    particle.pressure_force = Vec2::new(0.0, 0.0);

                    // if particle is along the edge of the window, apply a force to push it away
                    if row == 0 || row == particle_grid.num_rows - 1 || col == 0 || col == particle_grid.num_cols - 1 {
                        particle.pressure_force += calculate_boundary_force(&transform);
                    }

                    for other_entity in neighbours.iter() {
                         // We need to wrap this in "unsafe" because of the second "get" called on particles here.
                        let (other_particle, _, _) = particles.get_unchecked(*other_entity).unwrap();

                        // direction and distance from particle to other particle
                        let dir: Vec2;
                        let dist = Vec2::new(other_particle.predicted_position.x - particle.predicted_position.x, other_particle.predicted_position.y - particle.predicted_position.y).length();

                        if dist >= SMOOTHING_RADIUS {
                            continue;
                        }

                        if dist == 0.0 {
                            // is there a better way of doing this?
                            dir = get_random_direction();
                        } else {
                            dir = Vec2::new(other_particle.predicted_position.x - particle.predicted_position.x, other_particle.predicted_position.y - particle.predicted_position.y) / dist;
                        }

                        // calculate pressure force due to other particle
                        let pressure_force: Vec2 = calculate_force_between_two_particles(dir, &dist, particle.density, other_particle.density);
                        
                        // calcualte viscosity force
                        let viscosity_force: Vec2 = calculate_viscosity_between_two_particles(&dist, &particle, &other_particle);

                        particle.pressure_force += pressure_force + viscosity_force;
                    }
                }
            }
        }
    });
}

// Calculate the viscosity force acting on particle 1 due to particle 2.
fn calculate_viscosity_between_two_particles(dist: &f32, particle_1: &Particle, particle_2: &Particle) -> Vec2 {
    let influence = viscosity_smoothing_kernel(&SMOOTHING_RADIUS, &dist);
    return VISCOSITY_STRENGTH * (particle_2.velocity - particle_1.velocity) * influence;
}

// Returns force acting on particle 1 due to particle 2.
fn calculate_force_between_two_particles(dir: Vec2, dist: &f32, density_1: f32, density_2: f32) -> Vec2 {
    let slope = smoothing_kernel_derivative(&SMOOTHING_RADIUS, &dist);
    let shared_pressure_force = calculate_shared_pressure(&density_1, &density_2);

    if 0.0 != density_2 {
        return shared_pressure_force * dir * slope * PARTICLE_MASS / density_2;
    }

    Vec2::new(0.0, 0.0)
}

fn get_random_direction() -> Vec2 {
    let mut rng = rand::thread_rng();
    return Vec2::new(rng.gen(), rng.gen());
}

// Bounce off walls of window.
fn resolve_collisions(particle: &mut Particle, transform: &mut Transform) {
    let half_bound_size_x = WINDOW_WIDTH / 2.0 - PARTICLE_RADIUS;
    let half_bound_size_y = WINDOW_HEIGHT / 2.0 - PARTICLE_RADIUS;

    if transform.translation.x.abs() > half_bound_size_x {
        transform.translation.x = half_bound_size_x * sign(transform.translation.x);
        particle.velocity.x *= -COLLISION_DAMPING;
    }

    if transform.translation.y.abs() > half_bound_size_y {
        transform.translation.y = half_bound_size_y * sign(transform.translation.y);
        particle.velocity.y *= -COLLISION_DAMPING;
    }
}

// Returns the sign of a number.
fn sign(x: f32) -> f32 {
    if x > 0.0 {
        return 1.0;
    } else if x < 0.0 {
        return -1.0;
    } else {
        return 0.0;
    }
}

// Custom SPH function for viscosity force.
fn viscosity_smoothing_kernel(radius: &f32, dst: &f32) -> f32 {
    let value = f32::max(0., radius * radius - dst * dst);
    value * value * value / (radius * radius * radius * radius * radius * radius)
}

// Calculate the relative 'influence' of a particle at a given distance from a point.
fn smoothing_kernel(radius: &f32, dst: &f32) -> f32 {
    if dst >= radius {
        return 0.0;
    }

    let volume = (std::f32::consts::PI * radius.powi(4)) / 6.0;
    return (radius - dst) * (radius - dst) / volume;
}

// Calculates gradient of the smoothing kernel at a given distance.
fn smoothing_kernel_derivative(radius: &f32, dst: &f32) -> f32 {
    let scale = 12. / (std::f32::consts::PI * radius.powi(4));
    return scale * (dst - radius);
}

// Calculate and cache the density of each particle.
fn update_densities(mut particles: Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>, mut materials: ResMut<Assets<ColorMaterial>>, mut particle_grid: Query<&mut ParticleGrid>) {
    // Use partition grid to calculate density of each particle.
    let grid = particle_grid.single_mut();

    grid.particles.par_iter().enumerate().for_each(|(row, _)| {
        for col in 0..grid.num_cols {
            for entity in grid.particles[row][col].iter() {
                unsafe {
                    let (mut particle, _, _) = particles.get_unchecked(*entity).unwrap();
                    let neighbours = grid.get_particle_neighbours(&entity, row, col);
                    particle.density = 0.0;

                    for other_entity in neighbours.iter() {
                        // This line is the reason we need the "unsafe" above. 
                        // The borrow checker doesn't like the second "get" called on particles here.
                        let (other_particle, _, _) = particles.get_unchecked(*other_entity).unwrap();

                        let distance = Vec2::new(other_particle.predicted_position.x - particle.predicted_position.x, other_particle.predicted_position.y - particle.predicted_position.y).length();
                        let influence = smoothing_kernel(&SMOOTHING_RADIUS, &distance);
                        particle.density += PARTICLE_MASS * influence;
                    }
                }
            }
        }
    });

    // Set color based on velocity/density
    // Note: can't use par_iter_mut here because we need to access the materials asset,
    // and we can't access resources mutably without using locks.
    for (particle, _, mut color_mat) in particles.iter_mut() {
            let color_to_change: &mut Color;

            if let Some(color_material) = &mut color_mat.0 {
                color_to_change = &mut color_material.color;
            } else if let Some(color_handle) = color_mat.1 {
                let color_material = materials.get_mut(color_handle).unwrap();
                color_to_change = &mut color_material.color;
            } else {
                return;
            }

            set_particle_color_based_on_property(&particle, color_to_change);
        }
}

// red = density above target, green = at target density, blue = density below target
fn set_particle_color_based_on_property(particle: &Particle, color: &mut Color) {
    // how far from the target property (density/velocity) the particle is
    let dist_from_target_property: f32;
    let red: f32;
    let green: f32;

    if VISUALIZE_COLOR_BASED_ON == "density" {
        dist_from_target_property = (particle.density - TARGET_DENSITY) * 80.0;
    } else {
        // dist_from_target_property = (particle.velocity.length() * 0.7) - 0.5;
        dist_from_target_property = (particle.velocity.length() * 2.0) - 0.5;
    }

    if dist_from_target_property > 0.0 {
        red = f32::min(dist_from_target_property, 1.0);
        green = 0.0;
    } else {
        red = 0.0;
        green = f32::min(-dist_from_target_property, 1.0);
    }

    let blue: f32 = f32::min(1.0, 1. - (dist_from_target_property) - (red + green));
        
    *color = Color::rgb(red, green, blue);
}

fn calculate_shared_pressure(density_1: &f32, density_2: &f32) -> f32 {
    return (convert_density_to_pressure(density_1) + convert_density_to_pressure(density_2)) / 2.0;
}

// Converts a particle's density to a pressure force.
fn convert_density_to_pressure(density: &f32) -> f32 {
    return PRESSURE_MULTIPLIER * (density - TARGET_DENSITY);
}