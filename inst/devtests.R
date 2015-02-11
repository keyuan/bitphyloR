n0 <- Node$new()
tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1,
                 dpLambda = 1)
res1 <- tssb$GetMixture()
res2 <- tssb$ConvertTssbToIgraph()

plot(res2$g, layout = layout.reingold.tilford(res2$g),
     vertex.label = get.vertex.attribute(res2$g, name = "size"))

tt <- tssb$root

tssb$CullTree()

tt1 <- tssb$root
res3 <- tssb$ConvertTssbToIgraph()
plot(res3$g,
     layout = layout.reingold.tilford(res3$g),
     vertex.label = get.vertex.attribute(res3$g, name = "size"))
