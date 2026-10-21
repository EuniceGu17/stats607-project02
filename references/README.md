# Reference values

`table2.csv` transcribes the numerical entries of Table 2 (journal p. 486) of
Storey (2002), from https://genomics.princeton.edu/storeylab/papers/directfdr.pdf.
They are published, rounded simulation results, not exact theoretical values.
The transcription preserves the last BH power entry 0.133 as printed; do not
silently replace it with the current reproduction's result.

It also retains the printed 0.129 at pi0=0.8, gamma=0.001. That entry was
checked visually against the user's local journal PDF (journal p. 486),
and is a documented unresolved mismatch in the reproduction.

`paper_bh_threshold` is the mean BH rejection threshold, not its FDR level.
The original table does not report achieved BH FDR. No original code or seed
was used. Read the author's corrections alongside the paper:
https://genomics.princeton.edu/storeylab/papers/directfdr_corrections.pdf.
