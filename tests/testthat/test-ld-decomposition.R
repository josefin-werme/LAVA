# Equivalence of the "gram" and "svd" LD decompositions (decompose.ld, plink reference mode) and of their fallbacks.
# Downstream results depend on the decomposition R only through ncol(R) (K) and R %*% t(R), which are invariant to the
# arbitrary signs of the principal components, so those are the quantities compared.

decompose = function(G, method, prune.thresh=99, max.block.size=3000) {
	LAVA:::decompose.ld(list(mode="plink", data=G), prune.thresh, max.block.size, ld.decomposition=method)
}

# genotypes (0/1/2) for n individuals and p SNPs with LD decaying along the chromosome (AR(1) liabilities per haplotype)
simulate.genotypes = function(n, p, rho=0.9, seed=1) {
	set.seed(seed)
	haplotype = function() {
		z = matrix(rnorm(n * p), n, p)
		for (j in 2:p) z[,j] = rho * z[,j-1] + sqrt(1 - rho^2) * z[,j]
		maf = runif(p, 0.05, 0.5)
		1 * (pnorm(z) < rep(maf, each=n))
	}
	G = haplotype() + haplotype()
	storage.mode(G) = "double"
	G
}

expect_equivalent_decomposition = function(R.gram, R.svd, tol=1e-8) {
	expect_equal(ncol(R.gram), ncol(R.svd))
	P.gram = tcrossprod(R.gram); P.svd = tcrossprod(R.svd)
	expect_lt(norm(P.gram - P.svd, "F") / norm(P.svd, "F"), tol)
}

test_that("gram and svd decompositions are equivalent when the reference has more individuals than SNPs", {
	G = simulate.genotypes(n=1500, p=300)
	expect_equivalent_decomposition(decompose(G, "gram"), decompose(G, "svd"))
})

test_that("gram and svd decompositions are equivalent when the reference has fewer individuals than SNPs", {
	G = simulate.genotypes(n=250, p=800, seed=2)
	expect_equivalent_decomposition(decompose(G, "gram"), decompose(G, "svd"))
})

test_that("gram and svd decompositions are equivalent for a square genotype matrix", {
	G = simulate.genotypes(n=400, p=400, seed=3)
	expect_equivalent_decomposition(decompose(G, "gram"), decompose(G, "svd"))
})

test_that("the number of retained PCs agrees across pruning thresholds", {
	G1 = simulate.genotypes(n=1200, p=250, seed=4); G2 = simulate.genotypes(n=200, p=600, seed=5)
	for (thresh in c(50, 90, 95, 99, 99.9)) {
		expect_equal(ncol(decompose(G1, "gram", thresh)), ncol(decompose(G1, "svd", thresh)))
		expect_equal(ncol(decompose(G2, "gram", thresh)), ncol(decompose(G2, "svd", thresh)))
	}
})

test_that("missing genotypes are handled as in the svd route", {
	G = simulate.genotypes(n=800, p=200, seed=6)
	G[sample(length(G), 500)] = NA
	expect_equivalent_decomposition(decompose(G, "gram"), decompose(G, "svd"))
})

test_that("the Gram route declines (so decompose.ld uses svd) when a retained eigenvalue is below its tolerance", {
	G = simulate.genotypes(n=500, p=150, seed=7)
	X = scale(G); X[is.na(X)] = 0; X = scale(X)
	expect_false(is.null(LAVA:::decompose.gram(X, prune.thresh=99)))
	expect_null(LAVA:::decompose.gram(X, prune.thresh=99, rel.tol=1))	# lambda[K] > lambda[1] is impossible
})

test_that("gram falls back to svd (and its fallbacks) for non-finite standardised data", {
	G = simulate.genotypes(n=300, p=100, seed=8)
	G[,10] = 1	# monomorphic SNP: standardisation gives NaN
	X = scale(G); X[is.na(X)] = 0; X = scale(X)
	expect_null(LAVA:::decompose.gram(X, prune.thresh=99))
	run = function(method) tryCatch(decompose(G, method), error=function(e) conditionMessage(e))
	expect_identical(run("gram"), run("svd"))
})

test_that("pre-computed LD input (including block-wise decomposition) is unaffected by the method", {
	G = simulate.genotypes(n=600, p=240, seed=9)
	ld = list(mode="ld", ld=cor(G))
	for (block in c(3000, 100)) {
		expect_identical(LAVA:::decompose.ld(ld, 99, block, ld.decomposition="gram"),
		                 LAVA:::decompose.ld(ld, 99, block, ld.decomposition="svd"))
	}
})

test_that("process.locus gives equivalent results with both methods on the example data", {
	data.dir = test_path("..", "..", "vignettes", "data")
	skip_if_not(file.exists(file.path(data.dir, "g1000_test.bed")), "example data not available")
	info = read.table(file.path(data.dir, "input.info.txt"), header=T, stringsAsFactors=F)
	info$filename = file.path(data.dir, basename(info$filename))
	info.file = tempfile(fileext=".txt"); write.table(info, info.file, quote=F, row.names=F, sep="\t")
	input = suppressWarnings(process.input(info.file, file.path(data.dir, "sample.overlap.txt"), file.path(data.dir, "g1000_test"),
	                                       phenos=c("asthma", "diabetes", "hypothyroidism", "rheuma")))
	loci = read.loci(file.path(data.dir, "test.loci"))
	for (i in 1:2) {
		gram = process.locus(loci[i,], input)
		svd = process.locus(loci[i,], input, ld.decomposition="svd")
		expect_equal(gram$K, svd$K)
		for (field in c("h2.obs", "h2.latent", "sigma", "omega", "omega.cor")) expect_equal(gram[[field]], svd[[field]], tolerance=1e-8)
		expect_equal(crossprod(gram$delta), crossprod(svd$delta), tolerance=1e-8)
		expect_equal(run.univ(gram), run.univ(svd))
		set.seed(1); b.gram = run.bivar(gram, adap.thresh=NULL); set.seed(1); b.svd = run.bivar(svd, adap.thresh=NULL)
		expect_equal(b.gram, b.svd, tolerance=1e-6)
	}
	expect_error(process.locus(loci[1,], input, ld.decomposition="qr"))
})
