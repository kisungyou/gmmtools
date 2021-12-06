# # add ellipse -------------------------------------------------------------
# # https://cookierobotics.com/007/
# 
# # Generate 20 datasets with noise and fix GMM
# list_gmm  = list()
# for (i in 1:20){
#   data_i = T4cluster::genSMILEY(sd=0.25)$data
#   list_gmm[[i]] = gmm(data_i, k=12)
# }
# 
# # Find the average of models
# gcenter = wsum(list_gmm)
# 
# # Do Reduction using Two Methods
# grM = gmr2022Y(gcenter, merger="moment")
# 
# require("ggplot2")
# npts    = 200
# pgrid   = as.matrix(expand.grid(x=seq(from=-2,to=2,length.out=npts),
#                   y=seq(from=-2,to=2,length.out=npts)))
# probM   = gmmdensity(grM, data=pgrid)
# df_prob = data.frame(x=pgrid[,1], y=pgrid[,2], density=probM)
# 
# obj_plot = ggplot2::ggplot(df_prob, aes(x=x,y=y)) +  # ,z=density
#   geom_raster(aes(fill=density)) +
#   scale_fill_viridis_c() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   coord_fixed(xlim=c(-2,2), ylim=c(-2,2)) + 
#   theme(legend.position = "bottom")
# 
# # ellipse by for loop
# 
# vec_t = seq(from=0, to=2*pi, length.out=200)
# cos_t = cos(vec_t)
# sin_t = sin(vec_t)
# 
# K = length(grM$weight)
# for (k in 1:K){
#   tmp_mu  = as.vector(grM$mean[k,])
#   tmp_cov = as.matrix(grM$variance[,,k])
#   eig_cov = base::eigen(tmp_cov)
#   
#   # tmp_coord = t(eig_cov$vectors%*%rbind(sqrt(eig_cov$values[1])*cos_t, sqrt(eig_cov$values[2])*sin_t))
#   tmp_coord = t(eig_cov$vectors%*%rbind(sqrt(2.772589*eig_cov$values[1])*cos_t, sqrt(2.772589*eig_cov$values[2])*sin_t))
#   tmp_coord = sweep(tmp_coord, 2, tmp_mu, FUN="+")
#   tmp_df    = data.frame(xx=tmp_coord[,1], yy=tmp_coord[,2])
#   # tmp_df    = data.frame(xx=c(tmp_coord[,1], tmp_coord[1,1]), yy=c(tmp_coord[,2], tmp_coord[1,2]))
#   obj_plot  = obj_plot + 
#     geom_path(data=tmp_df, aes(x=xx, y=yy), color="white", size=1.5)
# }
# plot(obj_plot)