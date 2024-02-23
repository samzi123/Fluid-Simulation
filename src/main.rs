mod particle_grid;

use bevy::prelude::*;
use bevy::window::{PrimaryWindow, PresentMode};
use rand::Rng;
use particle_grid::{Particle, ParticleGrid};

/// We will store the world position of the mouse cursor here.
#[derive(Resource, Default)]
struct MyWorldCoords(Vec2);

/// Used to help identify our main camera
#[derive(Component)]
struct MainCamera;

const NUM_PARTICLES: usize = 4052;
const VISUALIZE_COLOR_BASED_ON: &str = "velocity"; // density or velocity
const PARTICLE_RADIUS: f32 = 4.;
const RESPOND_TO_MOUSE: bool = true;
const MOUSE_RADIUS: f32 = 150.0;
const GRAVITY: f32 = 2.0;
const WINDOW_WIDTH: f32 = 800.0;
const WINDOW_HEIGHT: f32 = 600.0;
const COLLISION_DAMPING: f32 = 0.6;
const PARTICLE_MASS: f32 = 1.0;
const SMOOTHING_RADIUS: f32 = 40.;
const TARGET_DENSITY: f32 = 0.0001;
// const PRESSURE_MULTIPLIER: f32 = 2000.0;
const PRESSURE_MULTIPLIER: f32 = 0.00000008;
const EPSILON: f32 = 1e-4;
const POLY6_CONSTANT: f32 = 315.0 / (64.0 * std::f32::consts::PI);

fn main() {
    App::new()
        .init_resource::<MyWorldCoords>()
        .insert_resource(ClearColor(Color::rgb(0.0, 0.0, 0.0)))
        .add_systems(Startup, setup)
        .add_systems(FixedUpdate, (update_particle_positions, update_densities))
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Fluid Simulation".into(),
                resolution: (WINDOW_WIDTH, WINDOW_HEIGHT).into(),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .run();
}

fn setup(
    mut commands: Commands,
    meshes: ResMut<Assets<Mesh>>,
    materials: ResMut<Assets<ColorMaterial>>,
    q_window: Query<&Window, With<PrimaryWindow>>,
) {
    let window_state = q_window.get_single().unwrap();
    commands.spawn((Camera2dBundle::default(), MainCamera));
    spawn_particles(commands, meshes, materials, window_state);
}

fn spawn_particles(
    mut commands: Commands,
    meshes: ResMut<Assets<Mesh>>,
    materials: ResMut<Assets<ColorMaterial>>,
    window_state: &Window,
) {
    // create spatial partitioning grid
    let mut particle_grid = ParticleGrid::new(SMOOTHING_RADIUS, (window_state.height() / SMOOTHING_RADIUS).ceil() as usize, (window_state.width() / SMOOTHING_RADIUS).ceil() as usize, window_state);
    particle_grid.spawn_particles(NUM_PARTICLES, PARTICLE_RADIUS, &mut commands, meshes, materials);
    commands.spawn(particle_grid);
}

// Updates the position of the particles based on their density and pressure.
fn update_particle_positions(
    mut particle_grid: Query<&mut ParticleGrid>,
    mut particles: Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>,
    mut mycoords: ResMut<MyWorldCoords>, 
    q_window: Query<&Window, With<PrimaryWindow>>, 
    q_camera: Query<(&Camera, &GlobalTransform), With<MainCamera>>, 
    time: Res<Time>,
) {
    let window_state = q_window.get_single().unwrap();
    let mut grid = particle_grid.single_mut();
    grid.update_particle_cells(&mut particles);
    
    // Reset velocity and pressure. If we don't do this, we get some weird momentum effects.
    for (mut particle, transform, _) in particles.iter_mut() {
        particle.predicted_position.x = f32::max(transform.translation.x + particle.velocity.x * time.delta_seconds(), -window_state.width() / 2.0 + PARTICLE_RADIUS);
        particle.predicted_position.x = f32::min(particle.predicted_position.x, window_state.width() / 2.0 - PARTICLE_RADIUS);
        particle.predicted_position.y = f32::max(transform.translation.y + particle.velocity.y * time.delta_seconds(), -window_state.height() / 2.0 + PARTICLE_RADIUS);
        particle.predicted_position.y = f32::min(particle.predicted_position.y, window_state.height() / 2.0 - PARTICLE_RADIUS);
        // particle.velocity = Vec2::new(0.0, 0.0);
    }

    if RESPOND_TO_MOUSE {
        let mouse_pos = get_mouse_world_position(&mut mycoords, &q_window, &q_camera);

        if mouse_pos != None {
            for (mut particle, transform, _) in particles.iter_mut() {
                let dir = mycoords.0 - transform.translation.xy();
                let euclidean_distance = dir.length();
                
                if euclidean_distance > MOUSE_RADIUS {
                    continue;
                }

                // particle.velocity -= smoothing_kernel(&MOUSE_RADIUS, &euclidean_distance) * dir * 40. * PRESSURE_MULTIPLIER * time.delta_seconds();
                particle.velocity -= smoothing_kernel(&MOUSE_RADIUS, &euclidean_distance) * dir * 10000. * time.delta_seconds();
            }
        }
    }

    update_pressure_forces(&mut particles, &grid);

    for (mut particle, mut transform, _) in particles.iter_mut() {
        // apply gravity
        // particle.velocity.y -= GRAVITY * time.delta_seconds();

        if particle.density != 0.0 {
            // F = m * a, so a = F / m
            let pressure_force = particle.pressure_force / particle.density;
            particle.velocity += pressure_force * time.delta_seconds();
        }
        
        transform.translation.x += particle.velocity.x;
        transform.translation.y += particle.velocity.y;

        resolve_collisions(&mut particle, &mut transform, window_state);
    }
}

// returns the mouse position in world coordinates
fn get_mouse_world_position(
    mycoords: &mut ResMut<MyWorldCoords>,
    // query to get the window (so we can read the current cursor position)
    q_window: &Query<&Window, With<PrimaryWindow>>,
    // query to get camera transform
    q_camera: &Query<(&Camera, &GlobalTransform), With<MainCamera>>,
) -> Option<Vec2> {
    // get the camera info and transform
    // assuming there is exactly one main camera entity, so Query::single() is OK
    let (camera, camera_transform) = q_camera.single();

    // There is only one primary window, so we can similarly get it from the query:
    let window = q_window.single();

    // check if the cursor is inside the window and get its position
    // then, ask bevy to convert into world coordinates, and truncate to discard Z
    if let Some(world_position) = window.cursor_position()
        .and_then(|cursor| camera.viewport_to_world(camera_transform, cursor))
        .map(|ray| ray.origin.truncate())
    {
        mycoords.0 = world_position;
        return Some(world_position);
    }

    // if the cursor is not inside the window, we don't update the world position
    return None;
}

// Calculate the force between all particles to simulate pressure.
fn update_pressure_forces(particles: &mut Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>, particle_grid: &ParticleGrid) {
    let cell_offsets: [(i32, i32); 9] = [
        (-1, -1),
        (-1, 0),
        (-1, 1),
        (0, -1),
        (0, 0),
        (0, 1),
        (1, -1),
        (1, 0),
        (1, 1),
    ];

    for row in 0..particle_grid.num_rows {
        for col in 0..particle_grid.num_cols {
            for entity in particle_grid.particles[row][col].iter() {
                unsafe {
                    let (mut particle, transform, _) = particles.get_unchecked(*entity).unwrap();
                    let neighbours = particle_grid.get_particle_neighbours(&entity, row, col);
                    particle.pressure_force = Vec2::new(0.0, 0.0);

                    for other_entity in neighbours.iter() {
                         // We need to wrap this in "unsafe" because of the second "get" called on particles here.
                        let (other_particle, other_transform, _) = particles.get_unchecked(*other_entity).unwrap();
                        let pressure_force = calculate_force_between_two_particles(particle.predicted_position, other_particle.predicted_position, particle.density, other_particle.density, &transform, &other_transform);
                        particle.pressure_force += pressure_force;
                    }
                }
            }
        }
    }
}

// Returns force acting on particle 1 due to particle 2.
fn calculate_force_between_two_particles(predicted_position_1: Vec2, predicted_position_2: Vec2, density_1: f32, density_2: f32, transform_1: &Transform, transform_2: &Transform) -> Vec2 {
    // let distance = Vec2::new(transform_2.translation.x - transform_1.translation.x, transform_2.translation.y - transform_1.translation.y).length();
    let distance = Vec2::new(predicted_position_2.x - predicted_position_1.x, predicted_position_2.y - predicted_position_1.y).length();
    let dir: Vec2;

    if distance > SMOOTHING_RADIUS {
        return Vec2::new(0.0, 0.0);
    } else if distance < EPSILON {
        // is there a better way of doing this?
        dir = get_random_direction();
    } else {
        // dir = Vec2::new(transform_2.translation.x - transform_1.translation.x, transform_2.translation.y - transform_1.translation.y) / distance;
        dir = Vec2::new(predicted_position_2.x - predicted_position_1.x, predicted_position_2.y - predicted_position_1.y) / distance;
    }

    let normalized_dir = dir.normalize();

    let slope = smoothing_kernel_derivative(&SMOOTHING_RADIUS, &distance);
    let shared_pressure_force = calculate_shared_pressure(&density_1, &density_2);

    if 0.0 != density_2 {
        // return shared_pressure_force * normalized_dir * slope * PARTICLE_MASS / density_2;
        // return shared_pressure_force * normalized_dir * slope * PARTICLE_MASS;
        let r_squared = distance * distance;
        let h_squared = SMOOTHING_RADIUS * SMOOTHING_RADIUS;
        
        if r_squared <= h_squared {
            let poly6_kernel = POLY6_CONSTANT * (h_squared - r_squared).powi(3);
            let pressure_force = poly6_kernel * slope * PARTICLE_MASS * PRESSURE_MULTIPLIER * (density_1 + density_2) / (2.0 * density_2);
            
            return pressure_force * normalized_dir;
        }
    }

    Vec2::new(0.0, 0.0)
}

fn get_random_direction() -> Vec2 {
    let mut rng = rand::thread_rng();
    return Vec2::new(rng.gen(), rng.gen());
}

// Bounce off walls of window.
fn resolve_collisions(particle: &mut Particle, transform: &mut Transform, window_state: &Window) {
    let half_bound_size_x = window_state.width() / 2.0 - PARTICLE_RADIUS;
    let half_bound_size_y = window_state.height() / 2.0 - PARTICLE_RADIUS;

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
    if dst >= radius {
        return 0.0;
    }

    let scale = 12. / (std::f32::consts::PI * radius.powi(4));
    return scale * (dst - radius);
}

// Calculate and cache the density of each particle.
fn update_densities(mut particles: Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>, mut materials: ResMut<Assets<ColorMaterial>>, mut particle_grid: Query<&mut ParticleGrid>) {
    // Use partition grid to calculate density of each particle.
    let grid = particle_grid.single_mut();
    let cell_offsets: [(i32, i32); 9] = [
        (-1, -1),
        (-1, 0),
        (-1, 1),
        (0, -1),
        (0, 0),
        (0, 1),
        (1, -1),
        (1, 0),
        (1, 1),
    ];

    for row in 0..grid.num_rows {
        for col in 0..grid.num_cols {
            for entity in grid.particles[row][col].iter() {
                unsafe {
                    let (mut particle, transform, _) = particles.get_unchecked(*entity).unwrap();
                    let neighbours = grid.get_particle_neighbours(&entity, row, col);
                    particle.density = 0.0;

                    for other_entity in neighbours.iter() {
                        // This line is the reason we need the "unsafe" above. 
                        // The borrow checker doesn't like the second "get" called on particles here.
                        let (_, other_transform, _) = particles.get_unchecked(*other_entity).unwrap();

                        let distance = Vec2::new(other_transform.translation.x - transform.translation.x, other_transform.translation.y - transform.translation.y).length();
                        let influence = smoothing_kernel(&SMOOTHING_RADIUS, &distance);
                        particle.density += PARTICLE_MASS * influence;
                    }
                }
            }
        }
    }

    // set color based on velocity/density
    for (particle, _, mut color_mat) in particles.iter_mut() {
        let color_to_change: &mut Color;

        if let Some(color_material) = &mut color_mat.0 {
            color_to_change = &mut color_material.color;
        } else if let Some(color_handle) = color_mat.1 {
            let color_material = materials.get_mut(color_handle).unwrap();
            color_to_change = &mut color_material.color;
        } else {
            continue;
        }

        set_particle_color_based_on_property(&particle, color_to_change);
    }
}

// red = density above target, green = at target density, blue = density below target
fn set_particle_color_based_on_property(particle: &Particle, color: &mut Color) {
    // takes some trial and error to get a value that looks good
    // let color_val = (particle.density) * 1700.;
    // how far from the target property (density/velocity) the particle is (positive or negative
    let dist_from_target_property: f32;

    if VISUALIZE_COLOR_BASED_ON == "density" {
        dist_from_target_property = (particle.density - TARGET_DENSITY) * 2000.0;
    } else {
        dist_from_target_property = (particle.velocity.length() * 0.7) - 0.5;
    }

    // let color_mutliplier: f32 = 2000.0;
    let red: f32;
    let blue: f32;

    if dist_from_target_property > 0.0 {
        red = f32::min(dist_from_target_property, 1.0);
        blue = 0.0;
    } else {
        red = 0.0;
        blue = f32::min(-dist_from_target_property, 1.0);
    }

    let green: f32 = f32::min(1.0, 1. - (dist_from_target_property) - (red + blue));
        
    let new_color = Color::rgb(red, green, blue);
    // let col = Color::rgb((particle_density * 100.), 0.0, 1. - (particle_density * 100.));
    *color = new_color;
}

fn calculate_shared_pressure(density_1: &f32, density_2: &f32) -> f32 {
    return (convert_density_to_pressure(density_1) + convert_density_to_pressure(density_2)) / 2.0;
}

// Converts a particle's density to a pressure force.
fn convert_density_to_pressure(density: &f32) -> f32 {
    return PRESSURE_MULTIPLIER * (density - TARGET_DENSITY);
}