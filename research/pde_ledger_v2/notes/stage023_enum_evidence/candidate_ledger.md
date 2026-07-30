# Candidate ledger — stage023 dimension-object enumeration

## Predicates declared before search results

- `P1 EXACT_TRIPLE_LITERAL`: every exact three-component numeric Wolfram list literal, with integer or rational components. This finds literal dimension vectors and deliberately accepts shape vectors and other false positives.
- `P2 DIM_LEXEME`: every identifier token containing case-insensitive `dim`. This finds dimension helpers, containers, local bindings, consumers, and deliberately accepts unrelated `Dimensions` calls.
- `P3 AUDIT_RECORD_KEY`: every exact string occurrence of one of the seven keys carried by `dimensionAudit`: `A0`, `A1`, `T0`, `T1`, `epsilon0`, `epsilon1`, `P0Physical`. This follows the record carriers without assuming that every occurrence is itself a dimension object.
- `P4 EXACT_NUMERIC_BINDING`: every identifier assignment to an exact integer or rational scalar. This finds candidate coefficients that may feed a walked expression and deliberately accepts counters, exponents, mutation locals, and print-text matches; disposition then asks whether the binding itself survives as an exponent-vector object.

All four predicates run over the complete 769-line source. Hits are individual regex matches, not deduplicated source lines. The locus includes a one-based column, and the matched text is the exact regex match.

Re-run from any working directory:

```bash
timeout 600 bash research/pde_ledger_v2/_scratch/stage023/enum/run_candidate_search.sh > research/pde_ledger_v2/_scratch/stage023/enum/candidate_ledger.md
```

## Raw search hits

```text
hit_id	predicate	locus	matched_text
P1-H001	P1	.wl:209:11	{0, 0, 0}
P1-H002	P1	.wl:238:8	{1, 0, 0}
P1-H003	P1	.wl:238:25	{1, 0, -1}
P1-H004	P1	.wl:238:46	{0, 0, -1}
P1-H005	P1	.wl:239:9	{0, 1, -1}
P1-H006	P1	.wl:239:27	{1, 1, -1}
P1-H007	P1	.wl:239:45	{0, 1, -1}
P1-H008	P1	.wl:239:63	{1, 1, -1}
P1-H009	P1	.wl:240:9	{-1, 1, -2}
P1-H010	P1	.wl:240:29	{0, 1, -2}
P1-H011	P1	.wl:240:49	{0, 1, -2}
P1-H012	P1	.wl:241:13	{0, 1, -2}
P1-H013	P1	.wl:241:34	{0, 1, -2}
P1-H014	P1	.wl:241:55	{0, 1, -2}
P1-H015	P1	.wl:242:13	{0, 0, -1}
P1-H016	P1	.wl:242:35	{0, 0, -1}
P1-H017	P1	.wl:242:55	{0, 0, -2}
P1-H018	P1	.wl:243:9	{-1/2, 1/2, -2}
P1-H019	P1	.wl:243:32	{-1/2, 1/2, -2}
P1-H020	P1	.wl:248:11	{0, 1, -1}
P1-H021	P1	.wl:248:31	{1, 1, -1}
P1-H022	P1	.wl:279:60	{1, 1, -1}
P1-H023	P1	.wl:283:63	{7, 0, 0}
P2-H001	P2	.wl:20:1	FAILDIMENSIONAL
P2-H002	P2	.wl:209:1	zeroDim
P2-H003	P2	.wl:210:1	dimScale
P2-H004	P2	.wl:212:1	dimOf
P2-H005	P2	.wl:212:54	dimensions
P2-H006	P2	.wl:214:41	zeroDim
P2-H007	P2	.wl:215:31	dims
P2-H008	P2	.wl:215:44	dims
P2-H009	P2	.wl:216:32	dimension
P2-H010	P2	.wl:217:33	dimOf
P2-H011	P2	.wl:217:42	dims
P2-H012	P2	.wl:221:47	dimension
P2-H013	P2	.wl:222:7	dimScale
P2-H014	P2	.wl:222:16	dimOf
P2-H015	P2	.wl:222:28	dims
P2-H016	P2	.wl:225:7	dimensions
P2-H017	P2	.wl:225:20	dimOf
P2-H018	P2	.wl:225:29	dims
P2-H019	P2	.wl:226:10	dimensions
P2-H020	P2	.wl:227:9	zeroDim
P2-H021	P2	.wl:228:36	dimensions
P2-H022	P2	.wl:229:17	dimension
P2-H023	P2	.wl:231:15	dimensions
P2-H024	P2	.wl:233:29	dimension
P2-H025	P2	.wl:237:1	baseDims
P2-H026	P2	.wl:244:14	zeroDim
P2-H027	P2	.wl:244:32	zeroDim
P2-H028	P2	.wl:244:50	zeroDim
P2-H029	P2	.wl:244:68	zeroDim
P2-H030	P2	.wl:247:1	expectedDims
P2-H031	P2	.wl:249:11	zeroDim
P2-H032	P2	.wl:249:28	zeroDim
P2-H033	P2	.wl:250:17	zeroDim
P2-H034	P2	.wl:250:40	zeroDim
P2-H035	P2	.wl:251:19	zeroDim
P2-H036	P2	.wl:254:1	dimensionAudit
P2-H037	P2	.wl:262:41	expectedDims
P2-H038	P2	.wl:262:61	dimOf
P2-H039	P2	.wl:262:86	dims
P2-H040	P2	.wl:266:57	expectedDims
P2-H041	P2	.wl:269:19	expectedDims
P2-H042	P2	.wl:270:6	DimensionalOk
P2-H043	P2	.wl:271:40	FAILDIMENSIONAL
P2-H044	P2	.wl:276:1	baselineDimAudit
P2-H045	P2	.wl:276:20	dimensionAudit
P2-H046	P2	.wl:277:3	baseDims
P2-H047	P2	.wl:279:1	corruptSourcedDims
P2-H048	P2	.wl:279:35	baseDims
P2-H049	P2	.wl:280:1	corruptSourcedDimAudit
P2-H050	P2	.wl:280:26	dimensionAudit
P2-H051	P2	.wl:281:3	corruptSourcedDims
P2-H052	P2	.wl:283:1	corruptFreeDims
P2-H053	P2	.wl:283:32	baseDims
P2-H054	P2	.wl:284:1	corruptFreeDimAudit
P2-H055	P2	.wl:284:23	dimensionAudit
P2-H056	P2	.wl:285:3	corruptFreeDims
P2-H057	P2	.wl:358:5	Dimensional
P2-H058	P2	.wl:358:19	FAILDIMENSIONAL
P2-H059	P2	.wl:381:54	dims
P2-H060	P2	.wl:381:60	dimAudit
P2-H061	P2	.wl:389:3	dims
P2-H062	P2	.wl:389:10	baseDims
P2-H063	P2	.wl:398:26	dims
P2-H064	P2	.wl:398:33	corruptSourcedDims
P2-H065	P2	.wl:415:3	dimAudit
P2-H066	P2	.wl:415:14	dimensionAudit
P2-H067	P2	.wl:416:5	dims
P2-H068	P2	.wl:431:6	Dimensional
P2-H069	P2	.wl:431:24	dimAudit
P2-H070	P2	.wl:431:34	DimensionalOk
P2-H071	P2	.wl:444:39	Dimension
P2-H072	P2	.wl:444:53	dimAudit
P2-H073	P2	.wl:509:41	dimensions
P2-H074	P2	.wl:509:55	Dimensions
P2-H075	P2	.wl:514:79	Dimensions
P2-H076	P2	.wl:607:1	runDimensionalGate
P2-H077	P2	.wl:607:33	dimFlag
P2-H078	P2	.wl:608:29	dimensional
P2-H079	P2	.wl:611:32	baselineDimAudit
P2-H080	P2	.wl:612:65	dimension
P2-H081	P2	.wl:612:77	baselineDimAudit
P2-H082	P2	.wl:614:5	expectedDims
P2-H083	P2	.wl:616:39	baselineDimAudit
P2-H084	P2	.wl:616:57	DimensionalOk
P2-H085	P2	.wl:617:39	corruptSourcedDimAudit
P2-H086	P2	.wl:618:44	corruptFreeDimAudit
P2-H087	P2	.wl:619:49	baselineDimAudit
P2-H088	P2	.wl:619:67	DimensionalOk
P2-H089	P2	.wl:620:82	corruptSourcedDimAudit
P2-H090	P2	.wl:620:117	FAILDIMENSIONAL
P2-H091	P2	.wl:621:34	dimension
P2-H092	P2	.wl:621:86	corruptFreeDimAudit
P2-H093	P2	.wl:622:58	corruptFreeDimAudit
P2-H094	P2	.wl:623:3	dimFlag
P2-H095	P2	.wl:623:54	dimension
P2-H096	P2	.wl:623:87	FAILDIMENSIONAL
P2-H097	P2	.wl:624:13	dimFlag
P2-H098	P2	.wl:682:84	Dimensions
P2-H099	P2	.wl:684:70	Dimensions
P2-H100	P2	.wl:686:21	dimensionAudit
P2-H101	P2	.wl:686:36	dims
P2-H102	P2	.wl:686:71	dimensional
P2-H103	P2	.wl:686:103	baselineDimAudit
P2-H104	P2	.wl:686:122	DimensionalOk
P2-H105	P2	.wl:693:5	baselineDimAudit
P2-H106	P2	.wl:693:23	corruptSourcedDimAudit
P2-H107	P2	.wl:693:47	corruptFreeDimAudit
P2-H108	P2	.wl:702:7	dimOf
P2-H109	P2	.wl:702:15	dimensionAudit
P2-H110	P2	.wl:739:36	dimensionalFlags
P2-H111	P2	.wl:745:3	dimensionalFlags
P2-H112	P2	.wl:745:22	runDimensionalGate
P2-H113	P2	.wl:748:34	dimensionalFlags
P3-H001	P3	.wl:116:5	"P0Physical"
P3-H002	P3	.wl:248:3	"A0"
P3-H003	P3	.wl:248:23	"A1"
P3-H004	P3	.wl:249:3	"T0"
P3-H005	P3	.wl:249:20	"T1"
P3-H006	P3	.wl:250:3	"epsilon0"
P3-H007	P3	.wl:250:26	"epsilon1"
P3-H008	P3	.wl:251:3	"P0Physical"
P3-H009	P3	.wl:257:5	"A0"
P3-H010	P3	.wl:257:17	"A1"
P3-H011	P3	.wl:257:29	"T0"
P3-H012	P3	.wl:257:41	"T1"
P3-H013	P3	.wl:258:5	"epsilon0"
P3-H014	P3	.wl:258:23	"epsilon1"
P3-H015	P3	.wl:258:41	"P0Physical"
P3-H016	P3	.wl:277:66	"P0Physical"
P3-H017	P3	.wl:281:76	"P0Physical"
P3-H018	P3	.wl:285:73	"P0Physical"
P3-H019	P3	.wl:416:48	"P0Physical"
P3-H020	P3	.wl:439:45	"T0"
P3-H021	P3	.wl:439:57	"T1"
P3-H022	P3	.wl:440:41	"A0"
P3-H023	P3	.wl:440:53	"A1"
P3-H024	P3	.wl:565:59	"T0"
P3-H025	P3	.wl:566:59	"T1"
P4-H001	P4	.wl:12:1	passCount = 0
P4-H002	P4	.wl:13:1	failCount = 0
P4-H003	P4	.wl:197:1	v0 = 1
P4-H004	P4	.wl:198:1	v1 = 1/2
P4-H005	P4	.wl:386:3	coeff1 = 1/2
P4-H006	P4	.wl:387:3	power1 = 3
P4-H007	P4	.wl:394:16	z0 = 0
P4-H008	P4	.wl:394:24	z1 = 0
P4-H009	P4	.wl:395:22	z0 = -2
P4-H010	P4	.wl:395:35	z1 = -2
P4-H011	P4	.wl:396:19	coeff1 = 1/3
P4-H012	P4	.wl:397:22	power1 = 2
P4-H013	P4	.wl:515:46	rank0=3
P4-H014	P4	.wl:565:30	T0=1/2
P4-H015	P4	.wl:566:30	T1=1/2
P4-H016	P4	.wl:675:51	ell=2
P4-H017	P4	.wl:719:39	v0=1
P4-H018	P4	.wl:719:44	v1=1/2
```

## Dispositions

```text
hit_id	disposition	row_key_or_code	reason
P1-H001	PROMOTE	zeroDim	source match supports this enumerated row
P1-H002	PROMOTE	sourced.a	source match supports this enumerated row
P1-H003	PROMOTE	sourced.cs	source match supports this enumerated row
P1-H004	PROMOTE	sourced.omega	source match supports this enumerated row
P1-H005	PROMOTE	sourced.M0	source match supports this enumerated row
P1-H006	PROMOTE	sourced.D1	source match supports this enumerated row
P1-H007	PROMOTE	sourced.R0	source match supports this enumerated row
P1-H008	PROMOTE	sourced.R1	source match supports this enumerated row
P1-H009	PROMOTE	sourced.D0	source match supports this enumerated row
P1-H010	PROMOTE	sourced.K0c	source match supports this enumerated row
P1-H011	PROMOTE	sourced.Keta	source match supports this enumerated row
P1-H012	PROMOTE	sourced.TOmega	source match supports this enumerated row
P1-H013	PROMOTE	sourced.Z0ret	source match supports this enumerated row
P1-H014	PROMOTE	sourced.Z1ret	source match supports this enumerated row
P1-H015	PROMOTE	sourced.OmegaU	source match supports this enumerated row
P1-H016	PROMOTE	sourced.OmegaW	source match supports this enumerated row
P1-H017	PROMOTE	sourced.Rmix	source match supports this enumerated row
P1-H018	PROMOTE	sourced.gU	source match supports this enumerated row
P1-H019	PROMOTE	sourced.gW	source match supports this enumerated row
P1-H020	PROMOTE	expected.A0	source match supports this enumerated row
P1-H021	PROMOTE	expected.A1	source match supports this enumerated row
P1-H022	PROMOTE	corruptSourcedDims	source match supports this enumerated row
P1-H023	PROMOTE	corruptFreeDims	source match supports this enumerated row
P2-H001	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H002	PROMOTE	zeroDim	source match supports this enumerated row
P2-H003	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H004	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H005	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H006	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H007	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H008	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H009	RULE_OUT	R9_OTHER	non-value lexical hit outside the membership rule
P2-H010	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H011	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H012	RULE_OUT	R9_OTHER	non-value lexical hit outside the membership rule
P2-H013	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H014	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H015	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H016	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H017	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H018	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H019	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H020	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H021	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H022	RULE_OUT	R9_OTHER	non-value lexical hit outside the membership rule
P2-H023	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H024	RULE_OUT	R9_OTHER	non-value lexical hit outside the membership rule
P2-H025	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H026	PROMOTE	sourced.etaNull	source match supports this enumerated row
P2-H027	PROMOTE	sourced.gain0	source match supports this enumerated row
P2-H028	PROMOTE	sourced.gain1	source match supports this enumerated row
P2-H029	PROMOTE	sourced.qfree	source match supports this enumerated row
P2-H030	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H031	PROMOTE	expected.T0	source match supports this enumerated row
P2-H032	PROMOTE	expected.T1	source match supports this enumerated row
P2-H033	PROMOTE	expected.epsilon0	source match supports this enumerated row
P2-H034	PROMOTE	expected.epsilon1	source match supports this enumerated row
P2-H035	PROMOTE	expected.P0Physical	source match supports this enumerated row
P2-H036	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H037	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H038	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H039	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H040	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H041	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H042	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H043	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H044	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H045	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H046	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H047	PROMOTE	corruptSourcedDims	source match supports this enumerated row
P2-H048	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H049	PROMOTE	corruptSourcedComputed	source match supports this enumerated row
P2-H050	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H051	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H052	PROMOTE	corruptFreeDims	source match supports this enumerated row
P2-H053	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H054	PROMOTE	corruptFreeComputed	source match supports this enumerated row
P2-H055	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H056	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H057	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H058	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H059	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H060	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H061	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H062	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H063	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H064	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H065	PROMOTE	caseComputedFamily	source match supports this enumerated row
P2-H066	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H067	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H068	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H069	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H070	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H071	RULE_OUT	R5_WRAPPER	association wrapper key whose contained computed set is counted
P2-H072	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H073	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H074	RULE_OUT	R4_SHAPE	Wolfram matrix-shape operator, not dimensional-analysis data
P2-H075	RULE_OUT	R4_SHAPE	Wolfram matrix-shape operator, not dimensional-analysis data
P2-H076	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H077	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H078	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H079	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H080	RULE_OUT	R9_OTHER	non-value lexical hit outside the membership rule
P2-H081	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H082	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H083	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H084	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H085	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H086	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H087	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H088	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H089	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H090	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H091	RULE_OUT	R9_OTHER	non-value lexical hit outside the membership rule
P2-H092	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H093	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H094	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H095	RULE_OUT	R9_OTHER	non-value lexical hit outside the membership rule
P2-H096	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H097	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H098	RULE_OUT	R4_SHAPE	Wolfram matrix-shape operator, not dimensional-analysis data
P2-H099	RULE_OUT	R4_SHAPE	Wolfram matrix-shape operator, not dimensional-analysis data
P2-H100	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H101	RULE_OUT	R2_LOCAL	parameter or local alias; the stable returned/bound object is counted instead
P2-H102	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H103	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H104	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H105	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H106	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H107	RULE_OUT	R6_REFERENCE	additional definition-container or consumer occurrence of an already counted canonical object
P2-H108	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H109	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H110	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H111	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P2-H112	RULE_OUT	R1_ROUTINE	function definition or call, not a dimension-valued binding
P2-H113	RULE_OUT	R3_STATUS	verdict token, Boolean status, or flag carrier rather than an exponent vector
P3-H001	RULE_OUT	R7_PHYSICAL	physical-expression key or physical-value access; its computed dimension is counted separately
P3-H002	PROMOTE	expected.A0	source match supports this enumerated row
P3-H003	PROMOTE	expected.A1	source match supports this enumerated row
P3-H004	PROMOTE	expected.T0	source match supports this enumerated row
P3-H005	PROMOTE	expected.T1	source match supports this enumerated row
P3-H006	PROMOTE	expected.epsilon0	source match supports this enumerated row
P3-H007	PROMOTE	expected.epsilon1	source match supports this enumerated row
P3-H008	PROMOTE	expected.P0Physical	source match supports this enumerated row
P3-H009	PROMOTE	baseline.A0	source match supports this enumerated row
P3-H010	PROMOTE	baseline.A1	source match supports this enumerated row
P3-H011	PROMOTE	baseline.T0	source match supports this enumerated row
P3-H012	PROMOTE	baseline.T1	source match supports this enumerated row
P3-H013	PROMOTE	baseline.epsilon0	source match supports this enumerated row
P3-H014	PROMOTE	baseline.epsilon1	source match supports this enumerated row
P3-H015	PROMOTE	baseline.P0Physical	source match supports this enumerated row
P3-H016	RULE_OUT	R6_REFERENCE	argument/consumer occurrence of an already counted computed-dimension object
P3-H017	RULE_OUT	R6_REFERENCE	argument/consumer occurrence of an already counted computed-dimension object
P3-H018	RULE_OUT	R6_REFERENCE	argument/consumer occurrence of an already counted computed-dimension object
P3-H019	RULE_OUT	R6_REFERENCE	argument/consumer occurrence of an already counted computed-dimension object
P3-H020	RULE_OUT	R7_PHYSICAL	physical-expression key or physical-value access; its computed dimension is counted separately
P3-H021	RULE_OUT	R7_PHYSICAL	physical-expression key or physical-value access; its computed dimension is counted separately
P3-H022	RULE_OUT	R7_PHYSICAL	physical-expression key or physical-value access; its computed dimension is counted separately
P3-H023	RULE_OUT	R7_PHYSICAL	physical-expression key or physical-value access; its computed dimension is counted separately
P3-H024	RULE_OUT	R7_PHYSICAL	physical-expression key or physical-value access; its computed dimension is counted separately
P3-H025	RULE_OUT	R7_PHYSICAL	physical-expression key or physical-value access; its computed dimension is counted separately
P4-H001	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H002	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H003	RULE_OUT	R8_TRANSIENT_NUMERIC	v0 evaluates to multiplicative identity and is absent from the walked A0 expression (.out:56)
P4-H004	RULE_OUT	R8_TRANSIENT_NUMERIC	v1 evaluates into the fused numeric atom I/2 in the walked A1 expression (.out:57)
P4-H005	RULE_OUT	R8_TRANSIENT_NUMERIC	coeff1 evaluates into a fused numeric atom before caseFor passes A1 to dimensionAudit (.wl:386,407,415-417)
P4-H006	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H007	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H008	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H009	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H010	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H011	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H012	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H013	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H014	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H015	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H016	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H017	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
P4-H018	RULE_OUT	R8_SCALAR_FALSE_POSITIVE	counter, mutation scalar, exponent, assertion text, or prose-string match; not a stable exponent-vector object
```

Rule-out codes: `R1_ROUTINE` function rather than value; `R2_LOCAL` local alias/parameter; `R3_STATUS` verdict or Boolean; `R4_SHAPE` matrix shape; `R5_WRAPPER` association wrapper; `R6_REFERENCE` repeated definition-container/consumer occurrence; `R7_PHYSICAL` physical value rather than dimension value; `R8_TRANSIENT_NUMERIC` scalar binding whose evaluated/fused numeric atom is transient rather than a stable exponent-vector object; `R8_SCALAR_FALSE_POSITIVE` counter/exponent/mutation/string false positive; `R9_OTHER` other non-value lexical hit.
