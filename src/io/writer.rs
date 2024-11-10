// use std::path::Path;
// use std::fs::File;
// use std::io::{self, Write, BufWriter};
// use anyhow::Result;

// pub struct FastxWriter {
//     writer: BufWriter<Box<dyn Write>>,
// }

// impl FastxWriter {
    
//     pub fn new<P: AsRef<Path>>(path: Option<P>) -> Result<Self> {
//         let writer = match path {
//             Some(path) => Box::new(File::create(path)?) as Box<dyn Write>,
//             None => Box::new(io::stdout()) as Box<dyn Write>,
//         };

//         Ok(FastxWriter {
//             writer: BufWriter::new(writer),
//         })
//     }

//     pub fn write_fasta(&mut self, id: &str, sequence: &[u8]) -> Result<()> {
//         writeln!(self.writer, ">{}", id)?;
//         writeln!(self.writer, "{}", String::from_utf8_lossy(sequence))?;
//         Ok(())
//     }

//     pub fn flush(&mut self) -> Result<()> {
//         self.writer.flush()?;
//         Ok(())
//     }
// }