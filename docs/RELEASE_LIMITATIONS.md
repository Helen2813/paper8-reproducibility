# Scientific and engineering scope of this candidate

1. **Leakage context:** the outcome-informed ComBat replay intentionally
   includes the 37 true test diagnoses. The derived classifier BA is not a
   prospectively unbiased clinical validation estimate.
2. **Data overlap:** identical CEL SHA-256 values establish identity of a
   recorded molecular measurement; they do not prove two GSM entries share
   a uniquely identified patient.
3. **Classifier convergence:** the historical frozen elastic-net parameters
   include `max_iter=100`. We reproduced this setting, and did not establish
   convergence to an optimized solution.
4. **Nested sample benchmark:** internal 324-probe cross-validation does not
   establish truly independent patient-level validation.
5. **Huang reconstruction:** five-gene random forest is a reconstruction
   inspired by the published approach, not the original fitted Huang model.
6. **Null models:** 18 historical conditions × 199 refits are described in
   the paper. The v3.70/v3.71 Python drivers were not recovered. Diagnostic
   pilots did not exactly match historical ComBat-null continuous scores.
   No software in this candidate purports to recreate all 3,582 original fits.
7. **Raw-data workflow:** standalone historical raw-CEL RMA/MAS5 checks passed
   on the author's workstation, but the candidate has no one-command
   path-independent GEO-to-all-results pipeline.
8. **Dependency reproducibility:** fresh environments passed on a single
   Apple Silicon Mac only. No second machine test or complete transitive
   package lock/guarantee of cross-platform reproducibility is included.
9. **Private inputs:** the saved-matrix core requires nine SHA-pinned files
   that are not distributed. This limits independent rerunnability.
10. **Permissions:** MIT covers only original code owned by the stated copyright holder. Manual upstream rights, coauthor code authorship,
   data sharing and person-level metadata disclosure review remain open.
