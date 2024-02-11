use bevy::{prelude::*, sprite::MaterialMesh2dBundle};
use bevy::window::{PrimaryWindow, PresentMode};
use rand::Rng;

/// We will store the world position of the mouse cursor here.
#[derive(Resource, Default)]
struct MyWorldCoords(Vec2);

/// Used to help identify our main camera
#[derive(Component)]
struct MainCamera;

#[derive(Component, Debug)]
struct Velocity {
    x: f32,
    y: f32,
}

#[derive(Component, Debug)]
struct Particle {
    velocity: Velocity,
    density: f32,
}

const NUM_PARTICLES: usize = 402;
const PARTICLE_RADIUS: f32 = 2.;
const RESPOND_TO_MOUSE: bool = false;
const GRAVITY: f32 = 9.81;
// const GRAVITY: f32 = 0.0;
const WINDOW_WIDTH: f32 = 800.0;
const WINDOW_HEIGHT: f32 = 600.0;
const COLLISION_DAMPING: f32 = 0.8;
const PARTICLE_MASS: f32 = 1.0;
const SMOOTHING_RADIUS: f32 = 40.;
const TARGET_DENSITY: f32 = 0.01;
const PRESSURE_MULTIPLIER: f32 = 100.0;

fn main() {
    App::new()
        .init_resource::<MyWorldCoords>()
        .insert_resource(ClearColor(Color::rgb(0.0, 0.0, 0.0)))
        .add_systems(Startup, setup)
        .add_systems(Update, update_position)
        .add_systems(Update, update_densities)
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
) {
    commands.spawn((Camera2dBundle::default(), MainCamera));
    spawn_particles(commands, meshes, materials);
}

fn spawn_particles(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    // spawn particles in a grid
    let particles_per_row: i32 = (NUM_PARTICLES as f32).sqrt() as i32;
    let particles_per_col: i32 = (NUM_PARTICLES as i32 - 1) / particles_per_row + 1;
    let spacing: f32 = PARTICLE_RADIUS * 6.0;

    for i in 0..NUM_PARTICLES {
        let x = ((i as i32 % particles_per_row - particles_per_row / 2) as f32 + 0.5) as f32 * spacing;
        let y = ((i as i32 / particles_per_row - particles_per_col / 2) as f32 + 0.5) as f32 * spacing;

        commands.spawn(
            Particle {
                velocity: Velocity { x: 0.0, y: 0.0 },
                density: 0.0,
            }
        ).insert(MaterialMesh2dBundle {
            mesh: meshes.add(shape::Circle::new(PARTICLE_RADIUS * 2.).into()).into(),
            material: materials.add(ColorMaterial::from(Color::WHITE)),
            transform: Transform::from_translation(Vec3::new(x, y, 0.)),
            ..default()
        });
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

// Adjust particle positions to even out their density.
fn update_position(
    mut particles: Query<(&mut Particle, &mut Transform)>, 
    mut mycoords: ResMut<MyWorldCoords>, 
    q_window: Query<&Window, With<PrimaryWindow>>, 
    q_camera: Query<(&Camera, &GlobalTransform), With<MainCamera>>, 
    time: Res<Time>,
) {
    let window_state = q_window.get_single().unwrap();

    if RESPOND_TO_MOUSE {
        let mouse_pos = get_mouse_world_position(&mut mycoords, &q_window, &q_camera);

        if mouse_pos == None {
            return;
        }
        
        for (mut particle, mut transform) in particles.iter_mut() {
            let distance = mycoords.0 - transform.translation.xy();
            let euclidean_distance = (distance.x * distance.x + distance.y * distance.y).sqrt();
            
            if euclidean_distance > 150.0 {
                continue;
            }

            particle.velocity.x = -distance.x * 0.1;
            particle.velocity.y = -distance.y * 0.1;

            // move towards/away from mouse
            transform.translation.x += particle.velocity.x * time.delta_seconds();
            transform.translation.y += particle.velocity.y * time.delta_seconds();
        }
    }

    // Reset velocities. If we don't do this, we get some weird momentum effects.
    for (mut particle, _) in particles.iter_mut() {
        particle.velocity = Velocity { x: 0.0, y: 0.0 };
    }

    update_pressure_forces(&mut particles, &time);

    for (mut particle, mut transform) in particles.iter_mut() {
        // apply gravity
        // particle.velocity.y -= GRAVITY * time.delta_seconds();

        transform.translation.x += particle.velocity.x;
        transform.translation.y += particle.velocity.y;

        resolve_collisions(&mut particle, &mut transform, window_state);
    }
}

fn update_pressure_forces(mut particles: &mut Query<(&mut Particle, &mut Transform)>, time: &Res<Time>) {
    let mut iter = particles.iter_combinations_mut();

    while let Some([(mut particle_1, mut transform_1), (mut particle_2, mut transform_2)]) = iter.fetch_next() {
        let distance = Vec2::new(transform_1.translation.x - transform_2.translation.x, transform_1.translation.y - transform_2.translation.y).length();
        let mut dir = Vec2::new(0.0, 0.0);

        // if distance > SMOOTHING_RADIUS {
        //     continue;
        // }
        
        if distance == 0.0 {
            dir = get_random_direction();
        } else {
            dir = Vec2::new(transform_1.translation.x - transform_2.translation.x, transform_1.translation.y - transform_2.translation.y) / distance;
        }

        let slope = smoothing_kernel_derivative(&SMOOTHING_RADIUS, &distance);
        
        let pressure_force = convert_density_to_pressure(&particle_1.density) * dir * slope * PARTICLE_MASS / particle_1.density;
        let pressure_acceleration = pressure_force / particle_1.density;
        particle_1.velocity.x += pressure_acceleration.x * time.delta_seconds();
        particle_1.velocity.y += pressure_acceleration.y * time.delta_seconds();
    }
}

fn get_random_direction() -> Vec2 {
    let mut rng = rand::thread_rng();
    let x: f32 = rng.gen();
    let y: f32 = rng.gen();
    return Vec2::new(x, y);
}

// Bounce off walls of window.
fn resolve_collisions(particle: &mut Particle, transform: &mut Transform, window_state: &Window) {
    let half_bound_size_x = window_state.width() / 2.0 - PARTICLE_RADIUS * 2.0;
    let half_bound_size_y = window_state.height() / 2.0 - PARTICLE_RADIUS * 2.0;

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

    // if dst > radius {
    //     return 0.0;
    // }

    // return f32::max(0.0, 1. - (dst / radius));
}

// Calculates gradient of the smoothing kernel at a given distance.
fn smoothing_kernel_derivative(radius: &f32, dst: &f32) -> f32 {
    if dst >= radius {
        return 0.0;
    }

    let scale = 12. / (std::f32::consts::PI * radius.powi(4));
    return scale * (dst - radius);
    // return *radius;
}

// Calculates the particle density of a given x and y coordinate.
fn calculate_density(sample_x: &f32, sample_y: &f32, particles: &Query<&Transform, With<Particle>>) -> f32 {
    let mut density: f32 = 0.0;

    for transform in particles.iter() {
        let distance = Vec2::new(sample_x - transform.translation.x, sample_y - transform.translation.y).length();
        let influence = smoothing_kernel(&SMOOTHING_RADIUS, &distance);
        density += PARTICLE_MASS * influence
    }

    density
}

// Calculate and cache the density of each particle.
fn update_densities(mut particles: Query<(&mut Particle, &Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>, transforms: Query<&Transform, With<Particle>>, mut materials: ResMut<Assets<ColorMaterial>> ) {
    for (mut particle, transform, color_mat) in particles.iter_mut() {
        particle.density = calculate_density(&transform.translation.x, &transform.translation.y, &transforms);
        // set color based on velocity
        let particle_velocity = particle.velocity.x.abs() + particle.velocity.y.abs();
        let particle_density = particle.density;
        let col = Color::rgb(particle_density / TARGET_DENSITY, 0.0, 1. - (particle_density / TARGET_DENSITY));

        if color_mat.0.is_some() {
            let some_color = &mut color_mat.0.unwrap().color;
            // *some_color = Color::rgb((particle_density * 100.), 0.0, 1. - (particle_density * 100.));
            *some_color = col;
        }
        if color_mat.1.is_some() {
            let some_handle = color_mat.1.unwrap();
            let some_color = &mut materials.get_mut(some_handle).unwrap().color;

            // *some_color = Color::rgb((particle_density * 100.), 0.0, 1. - (particle_density * 100.));
            *some_color = col;

        }
    }
}

fn calculate_pressure_force(
    sample_x: &f32,
    sample_y: &f32,
    particles: &Query<(&Particle, &Transform)>,
) -> Vec2 {
    let mut pressure_force = Vec2::new(0.0, 0.0);

    for (particle, transform) in particles.iter() {
        let distance = Vec2::new(sample_x - transform.translation.x, sample_y - transform.translation.y).length();
        
        if distance == 0.0 {
            continue;
        }

        let dir = Vec2::new(sample_x - transform.translation.x, sample_y - transform.translation.y) / distance;
        let slope = smoothing_kernel_derivative(&SMOOTHING_RADIUS, &distance);

        // gradient.x += influence * (sample_x - transform.translation.x);
        // gradient.y += influence * (sample_y - transform.translation.y);
        pressure_force += convert_density_to_pressure(&particle.density) * dir * slope * PARTICLE_MASS / particle.density;
    }

    return pressure_force;
}

fn convert_density_to_pressure(density: &f32) -> f32 {
    return PRESSURE_MULTIPLIER * (density - TARGET_DENSITY);
}