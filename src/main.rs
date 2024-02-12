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
struct Particle {
    velocity: Vec2,
    density: f32,
    pressure_force: Vec2,
}

const NUM_PARTICLES: usize = 402;
const PARTICLE_RADIUS: f32 = 4.;
const RESPOND_TO_MOUSE: bool = true;
const MOUSE_RADIUS: f32 = 150.0;
const GRAVITY: f32 = 50.1;
const WINDOW_WIDTH: f32 = 800.0;
const WINDOW_HEIGHT: f32 = 600.0;
const COLLISION_DAMPING: f32 = 0.8;
const PARTICLE_MASS: f32 = 1.0;
const SMOOTHING_RADIUS: f32 = 40.;
const TARGET_DENSITY: f32 = 0.0001;
const PRESSURE_MULTIPLIER: f32 = 1800.0;

fn main() {
    App::new()
        .init_resource::<MyWorldCoords>()
        .insert_resource(ClearColor(Color::rgb(0.0, 0.0, 0.0)))
        .add_systems(Startup, setup)
        .add_systems(Update, update_position)
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
                velocity: Vec2::new(0.0, 0.0),
                density: 0.0,
                pressure_force: Vec2::new(0.0, 0.0),
            }
        )
        // insert arrow to show velocity direction and circle for particle
        .insert(MaterialMesh2dBundle {
            mesh: meshes.add(shape::Circle {
                radius: PARTICLE_RADIUS,
                ..Default::default()
            }.into()).into(),
            material: materials.add(ColorMaterial::from(Color::WHITE)),
            transform: Transform::from_translation(Vec3::new(x, y, 0.0)),
            ..Default::default()
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
    mut particles: Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>,
    mut mycoords: ResMut<MyWorldCoords>, 
    q_window: Query<&Window, With<PrimaryWindow>>, 
    q_camera: Query<(&Camera, &GlobalTransform), With<MainCamera>>, 
    time: Res<Time>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    let window_state = q_window.get_single().unwrap();

    if RESPOND_TO_MOUSE {
        let mouse_pos = get_mouse_world_position(&mut mycoords, &q_window, &q_camera);

        if mouse_pos != None {
            for (mut particle, mut transform, _) in particles.iter_mut() {
                let dir = mycoords.0 - transform.translation.xy();
                let euclidean_distance = dir.length();
                
                if euclidean_distance > MOUSE_RADIUS {
                    continue;
                }

                particle.velocity -= smoothing_kernel(&MOUSE_RADIUS, &euclidean_distance) * dir * 60. * PRESSURE_MULTIPLIER;

                // move towards/away from mouse
                transform.translation.x += particle.velocity.x * time.delta_seconds();
                transform.translation.y += particle.velocity.y * time.delta_seconds();
            }
        }
    }

    // Reset velocity, densities, and pressure. If we don't do this, we get some weird momentum effects.
    for (mut particle, _, _) in particles.iter_mut() {
        particle.velocity = Vec2::new(0.0, 0.0);
        particle.pressure_force = Vec2::new(0.0, 0.0);
        particle.density = 0.0;
    }

    update_densities(&mut particles, &mut materials);
    update_pressure_forces(&mut particles);

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

// Calculate the force between all particles to simulate pressure.
fn update_pressure_forces(particles: &mut Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>) {
    let mut iter = particles.iter_combinations_mut();

    while let Some([(mut particle_1, transform_1, _), (mut particle_2, transform_2, _)]) = iter.fetch_next() {
        let distance = Vec2::new(transform_2.translation.x - transform_1.translation.x, transform_2.translation.y - transform_1.translation.y).length();
        let mut dir: Vec2;

        if distance == 0.0 {
            // is there a better way of doing this?
            dir = get_random_direction();
        } else {
            dir = Vec2::new(transform_2.translation.x - transform_1.translation.x, transform_2.translation.y - transform_1.translation.y) / distance;
        }

        let slope = smoothing_kernel_derivative(&SMOOTHING_RADIUS, &distance);
        
        if particle_2.density != 0.0 {
            let pressure_force = convert_density_to_pressure(&particle_2.density) * dir * slope * PARTICLE_MASS / particle_2.density;
            particle_1.pressure_force += pressure_force;
        }

        // do the same for particle 2 becuase iter_combinations_mut won't repeat this pair
        if particle_1.density != 0.0 {
            dir = -dir;
            let pressure_force_2 = convert_density_to_pressure(&particle_1.density) * dir * slope * PARTICLE_MASS / particle_1.density;
            particle_2.pressure_force += pressure_force_2;
        }
    }
}

fn get_random_direction() -> Vec2 {
    let mut rng = rand::thread_rng();
    return Vec2::new(rng.gen(), rng.gen());
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
fn update_densities(particles: &mut Query<(&mut Particle, &mut Transform, AnyOf<(&mut TextureAtlasSprite, &Handle<ColorMaterial>)>)>, materials: &mut ResMut<Assets<ColorMaterial>> ) {
    let mut iter = particles.iter_combinations_mut();

    while let Some([(mut particle_1, transform_1, _), (mut particle_2, transform_2, _)]) = iter.fetch_next() {
        let distance = Vec2::new(transform_2.translation.x - transform_1.translation.x, transform_2.translation.y - transform_1.translation.y).length();
        let influence = smoothing_kernel(&SMOOTHING_RADIUS, &distance);
        particle_1.density += PARTICLE_MASS * influence;
        particle_2.density += PARTICLE_MASS * influence;
    }

    // set color based on velocity/density
    for (particle, _, color_mat) in particles.iter_mut() {
        let particle_velocity = particle.velocity.x.abs() + particle.velocity.y.abs();
        // takes some trial and error to get a value that looks good
        let color_val = (particle.density) * 1700.;
        let col = Color::rgb(color_val, 0.0, 1. - color_val);
        // let col = Color::rgb((particle_density * 100.), 0.0, 1. - (particle_density * 100.));

        if color_mat.0.is_some() {
            let some_color = &mut color_mat.0.unwrap().color;
            *some_color = col;
        }
        if color_mat.1.is_some() {
            let some_handle = color_mat.1.unwrap();
            let some_color = &mut materials.get_mut(some_handle).unwrap().color;
            *some_color = col;
        }
    }
}

// Converts a particle's density to a pressure force.
fn convert_density_to_pressure(density: &f32) -> f32 {
    return PRESSURE_MULTIPLIER * (density - TARGET_DENSITY);
}