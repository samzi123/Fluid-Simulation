use bevy::{prelude::*, sprite::MaterialMesh2dBundle};
use bevy::window::{PrimaryWindow, PresentMode, WindowResized};

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

const NUM_PARTICLES: usize = 10;
const PARTICLE_RADIUS: f32 = 5.;
const RESPOND_TO_MOUSE: bool = false;
const GRAVITY: f32 = 9.81;
const WINDOW_WIDTH: f32 = 800.0;
const WINDOW_HEIGHT: f32 = 600.0;

fn main() {
    App::new()
        .init_resource::<MyWorldCoords>()
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
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    commands.spawn((Camera2dBundle::default(), MainCamera));

    for i in 0..NUM_PARTICLES {
        commands.spawn((
            MaterialMesh2dBundle {
                mesh: meshes.add(shape::Circle::new(PARTICLE_RADIUS * 2.).into()).into(),
                material: materials.add(ColorMaterial::from(Color::WHITE)),
                transform: Transform::from_translation(Vec3::new(-150. + (i * 5) as f32, 0. + (i * 5) as f32, 0.)),
                ..default()
            },
            Velocity { x: 0.0, y: 0.0 },
        ));
    }
}

// returns the mouse position in world coordinates
fn get_mouse_world_position(
    mycoords: &mut ResMut<MyWorldCoords>,
    // query to get the window (so we can read the current cursor position)
    q_window: Query<&Window, With<PrimaryWindow>>,
    // query to get camera transform
    q_camera: Query<(&Camera, &GlobalTransform), With<MainCamera>>,
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

// This system updates the position of the particles.
fn update_position(mut query: Query<(&mut Velocity, &mut Transform)>, mut mycoords: ResMut<MyWorldCoords>, q_window: Query<&Window, With<PrimaryWindow>>, q_camera: Query<(&Camera, &GlobalTransform), With<MainCamera>>, time: Res<Time>, window_query: Query<&Window, With<PrimaryWindow>>) {
    let window_state = window_query.get_single().unwrap();

    if RESPOND_TO_MOUSE {
        let mouse_pos = get_mouse_world_position(&mut mycoords, q_window, q_camera);

        if mouse_pos == None {
            return;
        }
        
        for (mut velocity, mut transform) in query.iter_mut() {
            let distance = mycoords.0 - transform.translation.xy();
            let euclidean_distance = (distance.x * distance.x + distance.y * distance.y).sqrt();
            
            if euclidean_distance > 150.0 {
                continue;
            }

            velocity.x = -distance.x * 0.1;
            velocity.y = -distance.y * 0.1;

            // move towards/away from mouse
            transform.translation.x += velocity.x * time.delta_seconds();
            transform.translation.y += velocity.y * time.delta_seconds();
        }
    } else {
        for (mut velocity, mut transform) in query.iter_mut() {
            // apply gravity
            velocity.y -= GRAVITY * time.delta_seconds();
            transform.translation.x += velocity.x;
            transform.translation.y += velocity.y;

            resolve_collisions(&mut velocity, &mut transform, window_state);
        }
    }
}

// Bounce off walls of window.
fn resolve_collisions(velocity: &mut Velocity, transform: &mut Transform, window_state: &Window) {
    let half_bound_size_x = window_state.width() / 2.0 - PARTICLE_RADIUS * 2.0;
    let half_bound_size_y = window_state.height() / 2.0 - PARTICLE_RADIUS * 2.0;

    if transform.translation.x.abs() > half_bound_size_x {
        transform.translation.x = half_bound_size_x * sign(transform.translation.x);
        velocity.x *= -1.0;
    }

    if transform.translation.y.abs() > half_bound_size_y {
        transform.translation.y = half_bound_size_y * sign(transform.translation.y);
        velocity.y *= -1.0;
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