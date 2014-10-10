This is a small parser and demonstration of the [GFA](http://lh3.github.io/2014/07/23/first-update-on-gfa/) format.

Usage:

	make
	./gfademo test.gfa

Expected output on test.gfa:

    ====== Link 1 + 2 + 5M =======
    M0: CGATGCAA.....
    M1: ...TGCAAAGTAC
    Matches: 5 Mismatches: 0 Gaps: 0 Identity: 1.000

    ====== Link 3 + 2 + 0M =======
    M0: TGCAACGTATAGACTTGTCAC..........
    M1: .....................TGCAAAGTAC
    Matches: 0 Mismatches: 0 Gaps: 0 Identity: 0.000

    ====== Link 3 + 4 - 1M1D2M1S =======
    M0: TGCAACGTATAGACTTGTCA....
    M1: ................G-CATATA
    Matches: 3 Mismatches: 0 Gaps: 1 Identity: 0.750

    ====== Link 4 - 5 + 0M =======
    M0: GCATATA........
    M1: .......CGATGATA
    Matches: 0 Mismatches: 0 Gaps: 0 Identity: 0.000

    ====== Containment 5 + 6 + 2 4M =======
    M0: CGATGA
    M1: ..ATGA
    Matches: 4 Mismatches: 0 Gaps: 0 Identity: 1.000
	

