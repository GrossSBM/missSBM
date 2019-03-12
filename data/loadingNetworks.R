# library(blockmodels)
# library(igraph)
#
# # American Political Books
# # https://wiki.cs.umd.edu/cmsc734_09/index.php?title=Social_Network_of_Political_Books
# load("AmericanPolBooks.Rdata")
# plot(IGpolbook,vertex.color = get.vertex.attribute(IGpolbook)$partition)
# modpolbook = BM_bernoulli("SBM_sym",as.matrix(get.adjacency(IGpolbook)))
# modpolbook$estimate()
# kbest = which.max(modpolbook$ICL)
# table(apply(modpolbook$memberships[[kbest]]$Z,1,which.max),get.vertex.attribute(IGpolbook)$partition)
#
#
# # American College Footbal
# # http://www.seldomusedreserve.com/?page_id=8805
# load("AmericanCollegeFootbal.Rdata")
# plot(IGcollFootball,vertex.color = as.numeric(as.factor(get.vertex.attribute(IGcollFootball)$partition)))
# # binarize
# A = as.matrix(get.adjacency(IGcollFootball))
# A[A>0] = 1
# modcollFootball = BM_bernoulli("SBM_sym",A)
# modcollFootball$estimate()
# kbest = which.max(modcollFootball$ICL)
# table(apply(modcollFootball$memberships[[kbest]]$Z,1,which.max),get.vertex.attribute(IGcollFootball)$partition)
#
#
# # French political blog
# # https://rdrr.io/cran/sand/man/fblog.html
# load("FrenchPolBlog.Rdata")
# FrenchPolBlog <- as.matrix(get.adjacency(IGblogFRpol))
# rownames(FrenchPolBlog) <- NULL
# colnames(FrenchPolBlog) <- NULL
#
# save(FrenchPolBlog, Id, file = "FrenchPolBlog.rda")
# plot(IGblogFRpol,vertex.color = as.numeric(as.factor(get.vertex.attribute(IGblogFRpol)$PolParty)))
# modFrench = BM_bernoulli("SBM_sym",as.matrix(get.adjacency(IGblogFRpol)))
# modFrench$estimate()
# kbest = which.max(modFrench$ICL)
# table(apply(modFrench$memberships[[kbest]]$Z,1,which.max),get.vertex.attribute(IGblogFRpol)$PolParty)
#
#
#
# # David Copperfield
# # http://konect.uni-koblenz.de/networks/adjnoun_adjacency
# load("DavidCopperfield.Rdata")
# plot(IGDavidCopperfield)
# modDC = BM_bernoulli("SBM_sym",as.matrix(get.adjacency(IGDavidCopperfield)))
# modDC$estimate()
# kbest = which.max(modDC$ICL)
# table(apply(modDC$memberships[[kbest]]$Z,1,which.max))
#
#
#
# # Eu research email
# # https://snap.stanford.edu/data/email-Eu-core.html
# load("Euresearch.Rdata")
# plot(IGEunet,vertex.color = get.vertex.attribute(IGEunet)$membership)
# modEu = BM_bernoulli("SBM",as.matrix(get.adjacency(IGEunet)))
# #modEu$estimate()
# #save(modEu,file="EuresBM.Rdata")
# load("EuresBM.Rdata")
# modEu$estimate()
# kbest = which.max(modEu$ICL)
# table(apply(modEu$memberships[[kbest]]$Z,1,which.max),get.vertex.attribute(IGEunet)$membership)
#
#
#
# # food web grassland
# # http://cosinproject.eu/extra/data/foodwebs/WEB.html
# load("foodweb.Rdata")
# plot(IGfoodweb)
# modFW = BM_bernoulli("SBM",as.matrix(get.adjacency(IGfoodweb)))
# modFW$estimate()
# kbest = which.max(modFW$ICL)
# table(apply(modFW$memberships[[kbest]]$Z,1,which.max))
