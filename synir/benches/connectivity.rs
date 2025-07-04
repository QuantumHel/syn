use criterion::{black_box, criterion_group, Criterion};
use rand::prelude::IndexedRandom;
use rand::seq::SliceRandom;
use rand::{rng, Rng};
use syn::architecture::connectivity::Connectivity;
use syn::architecture::Architecture;

fn random_connected_connectivity(
    num_nodes: usize,
    extra_edges: usize,
    subset_length: usize,
) -> (Connectivity, Vec<usize>, usize) {
    let mut rng = rng();
    let mut edges = Vec::new();

    let mut nodes: Vec<usize> = (0..num_nodes).collect();
    nodes.shuffle(&mut rng);
    for i in 1..num_nodes {
        edges.push((nodes[i - 1], nodes[i]));
    }

    // Add extra random edges:
    while edges.len() < (num_nodes - 1) + extra_edges {
        let a = rng.random_range(0..num_nodes);
        let b = rng.random_range(0..num_nodes);
        if a != b && !edges.contains(&(a, b)) && !edges.contains(&(b, a)) {
            edges.push((a, b));
        }
    }

    let mut subset_nodes = nodes.clone();
    subset_nodes.shuffle(&mut rng);
    let subset: Vec<usize> = subset_nodes.into_iter().take(subset_length).collect();

    let random_element = *subset.choose(&mut rng).expect("Subset cannot be empty");

    (Connectivity::from_edges(&edges), subset, random_element)
}

fn get_cx_ladder_connectivity((connectivity, nodes, root): &(Connectivity, Vec<usize>, usize)) {
    let _ = connectivity.get_cx_ladder(nodes, root);
}

pub fn connectivity_bench(c: &mut Criterion) {
    let input = random_connected_connectivity(30, 15, 10);
    c.bench_function("get_cx_ladder_connectivity: 30, 15, 10", |b| {
        b.iter(|| get_cx_ladder_connectivity(black_box(&input)))
    });

    let input = random_connected_connectivity(100, 50, 50);
    c.bench_function("get_cx_ladder_connectivity: 100, 15, 10", |b| {
        b.iter(|| get_cx_ladder_connectivity(black_box(&input)))
    });
}

criterion_group!(connectivity_benchmark, connectivity_bench);
