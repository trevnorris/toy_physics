import S11D5Invariants.ReconstructionBlock0
import S11D5Invariants.ReconstructionBlock1
import S11D5Invariants.ReconstructionBlock2
import S11D5Invariants.ReconstructionBlock3
import S11D5Invariants.ReconstructionBlock4
import S11D5Invariants.ReconstructionBlock5
import S11D5Invariants.ReconstructionBlock6
import S11D5Invariants.ReconstructionBlock7
import S11D5Invariants.ReconstructionBlock8
import S11D5Invariants.ReconstructionBlock9
import S11D5Invariants.ReconstructionBlock10
import S11D5Invariants.ReconstructionBlock11
import S11D5Invariants.ReconstructionBlock12
import S11D5Invariants.ReconstructionBlock13
import S11D5Invariants.ReconstructionBlock14
import S11D5Invariants.ReconstructionBlock15
import S11D5Invariants.ReconstructionBlock16

/-! Assemble the complete coefficient reconstruction. The full-group sufficiency
proof is separate; the finite certificate is independent of native Q9. -/
namespace S11D5Invariants
noncomputable section
set_option maxRecDepth 4096
set_option maxHeartbeats 3200000

theorem invariant_polynomial {Q : Quad} (hQ : OInvariant Q) (c : Coefficients)
    (hc : ∀ G, Q G = polynomial c (coordinates G)) :
    ∀ G, Q G = (c 6/2) * G.trace ^ 2 +
      (c 29/2) * (G*G).trace +
      c 25 * (G*G.transpose).trace := by
  obtain ⟨h0,h1,h2,h3,h4,h5,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20⟩ := reconstructedBlock0 hQ c hc
  obtain ⟨h21,h22,h23,h24,h26,h27,h28,h30,h31,h32,h33,h34,h35,h36,h37,h38,h39,h40,h41,h42⟩ := reconstructedBlock1 hQ c hc
  obtain ⟨h43,h44,h45,h46,h47,h48,h49,h50,h51,h52,h53,h54,h55,h56,h57,h58,h59,h60,h61,h62⟩ := reconstructedBlock2 hQ c hc
  obtain ⟨h63,h64,h65,h66,h67,h68,h69,h70,h71,h72,h73,h74,h75,h76,h77,h78,h79,h80,h81,h82⟩ := reconstructedBlock3 hQ c hc
  obtain ⟨h83,h84,h85,h86,h87,h88,h89,h90,h91,h92,h93,h94,h95,h96,h97,h98,h99,h100,h101,h102⟩ := reconstructedBlock4 hQ c hc
  obtain ⟨h103,h104,h105,h106,h107,h108,h109,h110,h111,h112,h113,h114,h115,h116,h117,h118,h119,h120,h121,h122⟩ := reconstructedBlock5 hQ c hc
  obtain ⟨h123,h124,h125,h126,h127,h128,h129,h130,h131,h132,h133,h134,h135,h136,h137,h138,h139,h140,h141,h142⟩ := reconstructedBlock6 hQ c hc
  obtain ⟨h143,h144,h145,h146,h147,h148,h149,h150,h151,h152,h153,h154,h155,h156,h157,h158,h159,h160,h161,h162⟩ := reconstructedBlock7 hQ c hc
  obtain ⟨h163,h164,h165,h166,h167,h168,h169,h170,h171,h172,h173,h174,h175,h176,h177,h178,h179,h180,h181,h182⟩ := reconstructedBlock8 hQ c hc
  obtain ⟨h183,h184,h185,h186,h187,h188,h189,h190,h191,h192,h193,h194,h195,h196,h197,h198,h199,h200,h201,h202⟩ := reconstructedBlock9 hQ c hc
  obtain ⟨h203,h204,h205,h206,h207,h208,h209,h210,h211,h212,h213,h214,h215,h216,h217,h218,h219,h220,h221,h222⟩ := reconstructedBlock10 hQ c hc
  obtain ⟨h223,h224,h225,h226,h227,h228,h229,h230,h231,h232,h233,h234,h235,h236,h237,h238,h239,h240,h241,h242⟩ := reconstructedBlock11 hQ c hc
  obtain ⟨h243,h244,h245,h246,h247,h248,h249,h250,h251,h252,h253,h254,h255,h256,h257,h258,h259,h260,h261,h262⟩ := reconstructedBlock12 hQ c hc
  obtain ⟨h263,h264,h265,h266,h267,h268,h269,h270,h271,h272,h273,h274,h275,h276,h277,h278,h279,h280,h281,h282⟩ := reconstructedBlock13 hQ c hc
  obtain ⟨h283,h284,h285,h286,h287,h288,h289,h290,h291,h292,h293,h294,h295,h296,h297,h298,h299,h300,h301,h302⟩ := reconstructedBlock14 hQ c hc
  obtain ⟨h303,h304,h305,h306,h307,h308,h309,h310,h311,h312,h313,h314,h315,h316,h317,h318,h319,h320,h321,h322⟩ := reconstructedBlock15 hQ c hc
  obtain ⟨h323,h324⟩ := reconstructedBlock16 hQ c hc
  intro G
  rw [hc]
  change polynomial c ![G 0 0,G 0 1,G 0 2,G 0 3,G 0 4,G 1 0,G 1 1,G 1 2,G 1 3,G 1 4,G 2 0,G 2 1,G 2 2,G 2 3,G 2 4,G 3 0,G 3 1,G 3 2,G 3 3,G 3 4,G 4 0,G 4 1,G 4 2,G 4 3,G 4 4] = _
  rw [polynomial_vec]
  simp only [h0, h1, h2, h3, h4, h5, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24, h26, h27, h28, h30, h31, h32, h33, h34, h35, h36, h37, h38, h39, h40, h41, h42, h43, h44, h45, h46, h47, h48, h49, h50, h51, h52, h53, h54, h55, h56, h57, h58, h59, h60, h61, h62, h63, h64, h65, h66, h67, h68, h69, h70, h71, h72, h73, h74, h75, h76, h77, h78, h79, h80, h81, h82, h83, h84, h85, h86, h87, h88, h89, h90, h91, h92, h93, h94, h95, h96, h97, h98, h99, h100, h101, h102, h103, h104, h105, h106, h107, h108, h109, h110, h111, h112, h113, h114, h115, h116, h117, h118, h119, h120, h121, h122, h123, h124, h125, h126, h127, h128, h129, h130, h131, h132, h133, h134, h135, h136, h137, h138, h139, h140, h141, h142, h143, h144, h145, h146, h147, h148, h149, h150, h151, h152, h153, h154, h155, h156, h157, h158, h159, h160, h161, h162, h163, h164, h165, h166, h167, h168, h169, h170, h171, h172, h173, h174, h175, h176, h177, h178, h179, h180, h181, h182, h183, h184, h185, h186, h187, h188, h189, h190, h191, h192, h193, h194, h195, h196, h197, h198, h199, h200, h201, h202, h203, h204, h205, h206, h207, h208, h209, h210, h211, h212, h213, h214, h215, h216, h217, h218, h219, h220, h221, h222, h223, h224, h225, h226, h227, h228, h229, h230, h231, h232, h233, h234, h235, h236, h237, h238, h239, h240, h241, h242, h243, h244, h245, h246, h247, h248, h249, h250, h251, h252, h253, h254, h255, h256, h257, h258, h259, h260, h261, h262, h263, h264, h265, h266, h267, h268, h269, h270, h271, h272, h273, h274, h275, h276, h277, h278, h279, h280, h281, h282, h283, h284, h285, h286, h287, h288, h289, h290, h291, h292, h293, h294, h295, h296, h297, h298, h299, h300, h301, h302, h303, h304, h305, h306, h307, h308, h309, h310, h311, h312, h313, h314, h315, h316, h317, h318, h319, h320, h321, h322, h323, h324]
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, Matrix.transpose_apply, sum_five]
  ring

end
end S11D5Invariants
