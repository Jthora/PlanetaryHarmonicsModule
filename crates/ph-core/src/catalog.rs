//! Event catalogues.
//!
//! Deliberately minimal: an event is a time, a location, and a size. Ingestion
//! for specific catalogues (Apollo PSE, USGS ComCat, tremor) builds on this.

/// A seismic event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Event {
    /// Ephemeris time, seconds past J2000 TDB.
    pub et: f64,
    /// Body-fixed latitude, degrees. `None` where a catalogue lacks locations —
    /// many Apollo deep moonquakes are identified only by nest.
    pub lat_deg: Option<f64>,
    /// Body-fixed longitude, degrees.
    pub lon_deg: Option<f64>,
    /// Depth in km, positive downward.
    pub depth_km: Option<f64>,
    /// Magnitude, on whatever scale the source catalogue uses.
    pub magnitude: Option<f64>,
}

/// A named set of events from one source.
#[derive(Debug, Clone, Default)]
pub struct Catalog {
    pub name: String,
    pub events: Vec<Event>,
}

impl Catalog {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            events: Vec::new(),
        }
    }

    /// Event times, ascending.
    pub fn times(&self) -> Vec<f64> {
        let mut t: Vec<f64> = self.events.iter().map(|e| e.et).collect();
        t.sort_by(|a, b| a.partial_cmp(b).unwrap());
        t
    }

    /// Events within `[start, end)`.
    pub fn window(&self, start: f64, end: f64) -> Catalog {
        Catalog {
            name: format!("{} [{start}, {end})", self.name),
            events: self
                .events
                .iter()
                .copied()
                .filter(|e| e.et >= start && e.et < end)
                .collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.events.len()
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(et: f64) -> Event {
        Event {
            et,
            lat_deg: None,
            lon_deg: None,
            depth_km: None,
            magnitude: None,
        }
    }

    #[test]
    fn times_come_back_sorted() {
        let mut c = Catalog::new("t");
        c.events = vec![ev(30.0), ev(10.0), ev(20.0)];
        assert_eq!(c.times(), vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn window_is_half_open() {
        let mut c = Catalog::new("t");
        c.events = vec![ev(10.0), ev(20.0), ev(30.0)];
        let w = c.window(10.0, 30.0);
        assert_eq!(w.len(), 2);
    }
}
