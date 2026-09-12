import S10Audit.CAS.PYLoci
import S10Audit.CAS.WLLoci

namespace S10Audit.CAS
open S10Pilot

theorem staticRank_locus_cross_engine (k : Vec 3) : PYLoci.staticRank k ↔ WLLoci.staticRank k := by
  rw [PYLoci.staticRank_geometry, WLLoci.staticRank_geometry]

theorem staticTransverse_locus_cross_engine (k : Vec 3) : PYLoci.staticTransverse k ↔ WLLoci.staticTransverse k := by
  rw [PYLoci.staticTransverse_geometry, WLLoci.staticTransverse_geometry]

theorem ordinaryRank_locus_cross_engine (k : Vec 3) : PYLoci.ordinaryRank k ↔ WLLoci.ordinaryRank k := by
  rw [PYLoci.ordinaryRank_geometry, WLLoci.ordinaryRank_geometry]

theorem ordinaryTransverse_locus_cross_engine (k : Vec 3) : PYLoci.ordinaryTransverse k ↔ WLLoci.ordinaryTransverse k := by
  rw [PYLoci.ordinaryTransverse_geometry, WLLoci.ordinaryTransverse_geometry]

theorem extraRank_locus_cross_engine (k : Vec 3) : PYLoci.extraRank k ↔ WLLoci.extraRank k := by
  rw [PYLoci.extraRank_geometry, WLLoci.extraRank_geometry]

theorem extraTransverse_locus_cross_engine (k : Vec 3) : PYLoci.extraTransverse k ↔ WLLoci.extraTransverse k := by
  rw [PYLoci.extraTransverse_geometry, WLLoci.extraTransverse_geometry]

end S10Audit.CAS
