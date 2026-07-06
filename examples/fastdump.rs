use astroimage::ImageAnalyzer;

fn main() {
    let file = std::env::args().nth(1).unwrap_or("tests/cocoon.fits".into());
    let cap: usize = std::env::args().nth(2).map(|v| v.parse().unwrap()).unwrap_or(500);
    let a = ImageAnalyzer::new().with_detection_sigma(5.0).with_max_stars(cap);
    let fast = a.detect_fast(&file).unwrap();
    let b = ImageAnalyzer::new().with_detection_sigma(5.0).with_max_stars(500);
    let slow = b.analyze(&file).unwrap();
    let fx: Vec<f32> = fast.stars.iter().map(|s| s.x).collect();
    let fy: Vec<f32> = fast.stars.iter().map(|s| s.y).collect();
    println!("fast: dims {}x{}  n={}  x[{:.0},{:.0}] y[{:.0},{:.0}]",
        fast.width, fast.height, fast.stars.len(),
        fx.iter().cloned().fold(f32::MAX,f32::min), fx.iter().cloned().fold(f32::MIN,f32::max),
        fy.iter().cloned().fold(f32::MAX,f32::min), fy.iter().cloned().fold(f32::MIN,f32::max));
    println!("slow: dims {}x{}  n={}", slow.width, slow.height, slow.stars.len());
    println!("--- top 8 slow (x y flux):");
    for s in slow.stars.iter().take(8) { println!("  {:8.2} {:8.2} {:12.0}", s.x, s.y, s.flux); }
    println!("--- top 8 fast (x y flux snr):");
    for s in fast.stars.iter().take(8) { println!("  {:8.2} {:8.2} {:12.0} {:8.1}", s.x, s.y, s.flux, s.snr); }
    // where do slow-top-100 stars sit in the fast list (rank by flux)?
    let mut ranks = Vec::new();
    for s in slow.stars.iter().take(100) {
        let mut best = f32::MAX; let mut bi = None;
        for (i, f) in fast.stars.iter().enumerate() {
            let d2 = (f.x-s.x).powi(2) + (f.y-s.y).powi(2);
            if d2 < best { best = d2; bi = Some(i); }
        }
        if best.sqrt() < 2.0 { ranks.push(bi.unwrap()); }
    }
    ranks.sort();
    println!("slow-top-100 matched in fast: {} ; fast-rank p10/p50/p90: {:?} {:?} {:?}",
        ranks.len(),
        ranks.get(ranks.len()/10), ranks.get(ranks.len()/2), ranks.get(ranks.len()*9/10));
    // nearest-fast distance for top 10 slow, and also with fast coords doubled
    let mut ds = Vec::new();
    for s in slow.stars.iter().take(100) {
        let mut best = f32::MAX;
        for f in fast.stars.iter() {
            let dx = f.x - s.x; let dy = f.y - s.y;
            best = best.min(dx*dx+dy*dy);
        }
        ds.push(best.sqrt());
    }
    ds.sort_by(|a,b| a.partial_cmp(b).unwrap());
    let n = ds.len();
    println!("nearest-dist p10={:.2} p50={:.2} p90={:.2} | coverage <2px {} <3px {} <4px {} <6px {} (of {})",
        ds[n/10], ds[n/2], ds[n*9/10],
        ds.iter().filter(|d| **d<2.0).count(), ds.iter().filter(|d| **d<3.0).count(),
        ds.iter().filter(|d| **d<4.0).count(), ds.iter().filter(|d| **d<6.0).count(), n);
}
