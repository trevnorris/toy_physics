# Phase A Existing-Provenance Cross-Check

You are running one blind Phase A modality. Do not read the outputs of other modalities.

Inputs:
- Stage list: `001, 002, 003, 004, 005, 006, 007, 008, 009, 010, 011, 012, 013, 014, 015, 016, 017, 018, 019, 020, 021, 022, 023, 024, 025, 026, 027, 028, 029, 030, 031, 032, 033, 034, 035, 036, 037, 038, 039, 040, 041, 042, 043, 044, 045, 046, 047, 048, 049, 050, 051, 052, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 082, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 093, 094, 095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253`
- Pass-2 reconciliation seed: `redteam/pass2/RECONCILIATION_AUGMENTATION.md`
- Pass-2 per-stage report paths: `- redteam/pass2/reports/stage_001.md
- redteam/pass2/reports/stage_002.md
- redteam/pass2/reports/stage_003.md
- redteam/pass2/reports/stage_004.md
- redteam/pass2/reports/stage_005.md
- redteam/pass2/reports/stage_006.md
- redteam/pass2/reports/stage_007.md
- redteam/pass2/reports/stage_008.md
- redteam/pass2/reports/stage_009.md
- redteam/pass2/reports/stage_010.md
- redteam/pass2/reports/stage_011.md
- redteam/pass2/reports/stage_012.md
- redteam/pass2/reports/stage_013.md
- redteam/pass2/reports/stage_014.md
- redteam/pass2/reports/stage_015.md
- redteam/pass2/reports/stage_016.md
- redteam/pass2/reports/stage_017.md
- redteam/pass2/reports/stage_018.md
- redteam/pass2/reports/stage_019.md
- redteam/pass2/reports/stage_020.md
- redteam/pass2/reports/stage_021.md
- redteam/pass2/reports/stage_022.md
- redteam/pass2/reports/stage_023.md
- redteam/pass2/reports/stage_024.md
- redteam/pass2/reports/stage_025.md
- redteam/pass2/reports/stage_026.md
- redteam/pass2/reports/stage_027.md
- redteam/pass2/reports/stage_028.md
- redteam/pass2/reports/stage_029.md
- redteam/pass2/reports/stage_030.md
- redteam/pass2/reports/stage_031.md
- redteam/pass2/reports/stage_032.md
- redteam/pass2/reports/stage_033.md
- redteam/pass2/reports/stage_034.md
- redteam/pass2/reports/stage_035.md
- redteam/pass2/reports/stage_036.md
- redteam/pass2/reports/stage_037.md
- redteam/pass2/reports/stage_038.md
- redteam/pass2/reports/stage_039.md
- redteam/pass2/reports/stage_040.md
- redteam/pass2/reports/stage_041.md
- redteam/pass2/reports/stage_042.md
- redteam/pass2/reports/stage_043.md
- redteam/pass2/reports/stage_044.md
- redteam/pass2/reports/stage_045.md
- redteam/pass2/reports/stage_046.md
- redteam/pass2/reports/stage_047.md
- redteam/pass2/reports/stage_048.md
- redteam/pass2/reports/stage_049.md
- redteam/pass2/reports/stage_050.md
- redteam/pass2/reports/stage_051.md
- redteam/pass2/reports/stage_052.md
- redteam/pass2/reports/stage_053.md
- redteam/pass2/reports/stage_054.md
- redteam/pass2/reports/stage_055.md
- redteam/pass2/reports/stage_056.md
- redteam/pass2/reports/stage_057.md
- redteam/pass2/reports/stage_058.md
- redteam/pass2/reports/stage_059.md
- redteam/pass2/reports/stage_060.md
- redteam/pass2/reports/stage_061.md
- redteam/pass2/reports/stage_062.md
- redteam/pass2/reports/stage_063.md
- redteam/pass2/reports/stage_064.md
- redteam/pass2/reports/stage_065.md
- redteam/pass2/reports/stage_066.md
- redteam/pass2/reports/stage_067.md
- redteam/pass2/reports/stage_068.md
- redteam/pass2/reports/stage_069.md
- redteam/pass2/reports/stage_070.md
- redteam/pass2/reports/stage_071.md
- redteam/pass2/reports/stage_072.md
- redteam/pass2/reports/stage_073.md
- redteam/pass2/reports/stage_074.md
- redteam/pass2/reports/stage_075.md
- redteam/pass2/reports/stage_076.md
- redteam/pass2/reports/stage_077.md
- redteam/pass2/reports/stage_078.md
- redteam/pass2/reports/stage_079.md
- redteam/pass2/reports/stage_080.md
- redteam/pass2/reports/stage_081.md
- redteam/pass2/reports/stage_082.md
- redteam/pass2/reports/stage_083.md
- redteam/pass2/reports/stage_084.md
- redteam/pass2/reports/stage_085.md
- redteam/pass2/reports/stage_086.md
- redteam/pass2/reports/stage_087.md
- redteam/pass2/reports/stage_088.md
- redteam/pass2/reports/stage_089.md
- redteam/pass2/reports/stage_090.md
- redteam/pass2/reports/stage_091.md
- redteam/pass2/reports/stage_092.md
- redteam/pass2/reports/stage_093.md
- redteam/pass2/reports/stage_094.md
- redteam/pass2/reports/stage_095.md
- redteam/pass2/reports/stage_096.md
- redteam/pass2/reports/stage_097.md
- redteam/pass2/reports/stage_098.md
- redteam/pass2/reports/stage_099.md
- redteam/pass2/reports/stage_100.md
- redteam/pass2/reports/stage_101.md
- redteam/pass2/reports/stage_102.md
- redteam/pass2/reports/stage_103.md
- redteam/pass2/reports/stage_104.md
- redteam/pass2/reports/stage_105.md
- redteam/pass2/reports/stage_106.md
- redteam/pass2/reports/stage_107.md
- redteam/pass2/reports/stage_108.md
- redteam/pass2/reports/stage_109.md
- redteam/pass2/reports/stage_110.md
- redteam/pass2/reports/stage_111.md
- redteam/pass2/reports/stage_112.md
- redteam/pass2/reports/stage_113.md
- redteam/pass2/reports/stage_114.md
- redteam/pass2/reports/stage_115.md
- redteam/pass2/reports/stage_116.md
- redteam/pass2/reports/stage_117.md
- redteam/pass2/reports/stage_118.md
- redteam/pass2/reports/stage_119.md
- redteam/pass2/reports/stage_120.md
- redteam/pass2/reports/stage_121.md
- redteam/pass2/reports/stage_122.md
- redteam/pass2/reports/stage_123.md
- redteam/pass2/reports/stage_124.md
- redteam/pass2/reports/stage_125.md
- redteam/pass2/reports/stage_126.md
- redteam/pass2/reports/stage_127.md
- redteam/pass2/reports/stage_128.md
- redteam/pass2/reports/stage_129.md
- redteam/pass2/reports/stage_130.md
- redteam/pass2/reports/stage_131.md
- redteam/pass2/reports/stage_132.md
- redteam/pass2/reports/stage_133.md
- redteam/pass2/reports/stage_134.md
- redteam/pass2/reports/stage_135.md
- redteam/pass2/reports/stage_136.md
- redteam/pass2/reports/stage_137.md
- redteam/pass2/reports/stage_138.md
- redteam/pass2/reports/stage_139.md
- redteam/pass2/reports/stage_140.md
- redteam/pass2/reports/stage_141.md
- redteam/pass2/reports/stage_142.md
- redteam/pass2/reports/stage_143.md
- redteam/pass2/reports/stage_144.md
- redteam/pass2/reports/stage_145.md
- redteam/pass2/reports/stage_146.md
- redteam/pass2/reports/stage_147.md
- redteam/pass2/reports/stage_148.md
- redteam/pass2/reports/stage_149.md
- redteam/pass2/reports/stage_150.md
- redteam/pass2/reports/stage_151.md
- redteam/pass2/reports/stage_152.md
- redteam/pass2/reports/stage_153.md
- redteam/pass2/reports/stage_154.md
- redteam/pass2/reports/stage_155.md
- redteam/pass2/reports/stage_156.md
- redteam/pass2/reports/stage_157.md
- redteam/pass2/reports/stage_158.md
- redteam/pass2/reports/stage_159.md
- redteam/pass2/reports/stage_160.md
- redteam/pass2/reports/stage_161.md
- redteam/pass2/reports/stage_162.md
- redteam/pass2/reports/stage_163.md
- redteam/pass2/reports/stage_164.md
- redteam/pass2/reports/stage_165.md
- redteam/pass2/reports/stage_166.md
- redteam/pass2/reports/stage_167.md
- redteam/pass2/reports/stage_168.md
- redteam/pass2/reports/stage_169.md
- redteam/pass2/reports/stage_170.md
- redteam/pass2/reports/stage_171.md
- redteam/pass2/reports/stage_172.md
- redteam/pass2/reports/stage_173.md
- redteam/pass2/reports/stage_174.md
- redteam/pass2/reports/stage_175.md
- redteam/pass2/reports/stage_176.md
- redteam/pass2/reports/stage_177.md
- redteam/pass2/reports/stage_178.md
- redteam/pass2/reports/stage_179.md
- redteam/pass2/reports/stage_180.md
- redteam/pass2/reports/stage_181.md
- redteam/pass2/reports/stage_182.md
- redteam/pass2/reports/stage_183.md
- redteam/pass2/reports/stage_184.md
- redteam/pass2/reports/stage_185.md
- redteam/pass2/reports/stage_186.md
- redteam/pass2/reports/stage_187.md
- redteam/pass2/reports/stage_188.md
- redteam/pass2/reports/stage_189.md
- redteam/pass2/reports/stage_190.md
- redteam/pass2/reports/stage_191.md
- redteam/pass2/reports/stage_192.md
- redteam/pass2/reports/stage_193.md
- redteam/pass2/reports/stage_194.md
- redteam/pass2/reports/stage_195.md
- redteam/pass2/reports/stage_196.md
- redteam/pass2/reports/stage_197.md
- redteam/pass2/reports/stage_198.md
- redteam/pass2/reports/stage_199.md
- redteam/pass2/reports/stage_200.md
- redteam/pass2/reports/stage_201.md
- redteam/pass2/reports/stage_202.md
- redteam/pass2/reports/stage_203.md
- redteam/pass2/reports/stage_204.md
- redteam/pass2/reports/stage_205.md
- redteam/pass2/reports/stage_206.md
- redteam/pass2/reports/stage_207.md
- redteam/pass2/reports/stage_208.md
- redteam/pass2/reports/stage_209.md
- redteam/pass2/reports/stage_210.md
- redteam/pass2/reports/stage_211.md
- redteam/pass2/reports/stage_212.md
- redteam/pass2/reports/stage_213.md
- redteam/pass2/reports/stage_214.md
- redteam/pass2/reports/stage_215.md
- redteam/pass2/reports/stage_216.md
- redteam/pass2/reports/stage_217.md
- redteam/pass2/reports/stage_218.md
- redteam/pass2/reports/stage_219.md
- redteam/pass2/reports/stage_220.md
- redteam/pass2/reports/stage_221.md
- redteam/pass2/reports/stage_222.md
- redteam/pass2/reports/stage_223.md
- redteam/pass2/reports/stage_224.md
- redteam/pass2/reports/stage_225.md
- redteam/pass2/reports/stage_226.md
- redteam/pass2/reports/stage_227.md
- redteam/pass2/reports/stage_228.md
- redteam/pass2/reports/stage_229.md
- redteam/pass2/reports/stage_230.md
- redteam/pass2/reports/stage_231.md
- redteam/pass2/reports/stage_232.md
- redteam/pass2/reports/stage_233.md
- redteam/pass2/reports/stage_234.md
- redteam/pass2/reports/stage_235.md
- redteam/pass2/reports/stage_236.md
- redteam/pass2/reports/stage_237.md
- redteam/pass2/reports/stage_238.md
- redteam/pass2/reports/stage_239.md
- redteam/pass2/reports/stage_240.md
- redteam/pass2/reports/stage_241.md
- redteam/pass2/reports/stage_242.md
- redteam/pass2/reports/stage_243.md
- redteam/pass2/reports/stage_244.md
- redteam/pass2/reports/stage_245.md
- redteam/pass2/reports/stage_246.md
- redteam/pass2/reports/stage_247.md
- redteam/pass2/reports/stage_248.md
- redteam/pass2/reports/stage_249.md
- redteam/pass2/reports/stage_250.md
- redteam/pass2/reports/stage_251.md
- redteam/pass2/reports/stage_252.md
- redteam/pass2/reports/stage_253.md`
- Checkpoint provenance seed: `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`

Task:
Every pass-2 value-reconciliation entry and checkpoint-provenance value that maps to these stages must become either a fit-insertion-point candidate or an explicit pure-identity classification. Do not exempt a value because it was previously audited.

Emit only YAML:

```yaml
modality: existing_provenance
candidates:
  - candidate_key:
    anchor_stage:
    parameter_names: []
    citation:
      path:
      line:
      excerpt:
    reason:
pure_identities:
  - stage:
    value:
    citation:
      path:
      line:
      excerpt:
    reason:
```
