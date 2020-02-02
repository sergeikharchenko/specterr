folder.name <- "spec.data "
dir.create(file.path (folder.name), showWarnings = F)
wind.size = 100; crop.size = 500; mov.step = 5L; tr.deg = 1; prc = 0.05; NLR = T; DEM.path = paste0(file.path (folder.name),"/DEM.tif"); border.path = paste0(file.path (folder.name), "/border.shp")

#####
### LOAD LIBRARIES
library(raster)
library(rgdal)
library(spectral)
library(spatial)
library(foreach)
library(doParallel)

#####
### READ RASTER DEM
gt <- raster(DEM.path) # считываем исходный растр

##### CREATE WORKSPACE AND TEMP FOLDERS
###
dir.create(path = paste0(getwd(),"/",file.path (folder.name, wind.size)), showWarnings = F) # Создаем папку для результатов

#####
### SET CROPING AND MOVING WINDOW OPTIONS
pix.size1 <- res(gt)[1] # разрешение растра по X
pix.size2 <- res(gt)[2] # разрешение растра по Y

gt.stack <- stack(raster(matrix(0, crop.size, crop.size))) # растр-стек для нарезки
ext.list <- list()

i.tiles <- ceiling(gt@ncols / (crop.size - wind.size + mov.step)) # столбцов в нарезке
j.tiles <- ceiling(gt@nrows / (crop.size - wind.size + mov.step)) # строк в нарезке

### Создаем X-ы и Y-и границ тайлов
xxx <- matrix(NA, i.tiles, 2)
for (i in 1:i.tiles) {
  xxx[i,] <- c(extent(gt)[1] + (i - 1) * (crop.size - (wind.size - mov.step))*res(gt)[1], 
               extent(gt)[1] + (i - 1) * (crop.size - (wind.size - mov.step))*res(gt)[1] + res(gt)[1] * crop.size)
}
yyy <- matrix(NA, j.tiles,2)
for (i in 1:j.tiles) {
  yyy[i,] <- c(extent(gt)[4] - (i - 1) * (crop.size - (wind.size - mov.step))*res(gt)[1], 
               extent(gt)[4] - (i - 1) * (crop.size - (wind.size - mov.step))*res(gt)[1] - res(gt)[1] * crop.size)
}
xy <- matrix(NA, i.tiles*j.tiles, 4)
for (i in 1:i.tiles) {
  for (j in 1:j.tiles) {
    xy[j.tiles*(i-1) + j,] <- c(xxx[i,], yyy[j,2:1])
  }
}

# Преобразуем X-ы и Y-и в векторные полигоны
tt <- lapply(apply(X = xy, MARGIN = 1, FUN = function(x) extent(x)), FUN = function(x) as(x, "SpatialPolygons"))
tt <- do.call(bind, tt)
tt$ID <- as.vector(t(matrix(1:length(tt), j.tiles, i.tiles)))
crs(tt) <- crs(gt)
plot(gt) ; plot(tt, add = T)
# writeOGR(tt, "tt.shp", "tt.shp", "ESRI Shapefile")

# Нарезаем ЦМР по полигонам
print(paste0("croping raster into mosaic with the ",i.tiles*j.tiles," segments"))
for (i in 1:i.tiles) {
  for (j in 1:j.tiles){
    if (extent(tt[j.tiles*(i-1) + j,])[2] <= extent(gt)[2]) {
      if (extent(tt[j.tiles*(i-1) + j,])[3] >= extent(gt)[3]) {
        raster.temp <- crop(gt, tt[j.tiles*(i-1) + j,])
      } else {
        raster.temp <- crop(gt, tt[j.tiles*(i-1) + j,]) # tt$ID[j.tiles*(i-1) + j],
        mNA <- matrix(NA, crop.size, crop.size)
        mNA[1:nrow(raster.temp), 1:ncol(raster.temp)] <- as.matrix(raster.temp)
        raster.temp <- raster(mNA)
      }
    } else {
      raster.temp <- crop(gt, tt[j.tiles*(i-1) + j,]) # tt$ID[j.tiles*(i-1) + j],
      mNA <- matrix(NA, crop.size, crop.size)
      mNA[1:nrow(raster.temp), 1:ncol(raster.temp)] <- as.matrix(raster.temp)
      raster.temp <- raster(mNA)
    }
    # plot(raster.temp, main = j.tiles*(i-1) + j)
    extent(raster.temp) <- extent(c(0,1,0,1))
    gt.stack[[j.tiles*(i-1) + j]] <- raster.temp
    print(j.tiles*(i-1) + j)
  }
}

# Создаем переменную с итоговыми экстентами моделей СХР, полученных из каждого растра ЦМР
xy2 <- xy
xy2[,c(1,3)] <- xy2[,c(1,3)] + res(gt)[1]*(wind.size - mov.step) / 2
xy2[,c(2,4)] <- xy2[,c(2,4)] - res(gt)[1]*(wind.size - mov.step) / 2
tt.ext <- apply(X = xy2, MARGIN = 1, FUN = function(x) extent(x))

#####
### CREATE INDICIES VECTORS
tiles <- i.tiles*j.tiles
ind.st <- seq(1, crop.size - wind.size, mov.step) + floor((crop.size - max(seq(1, crop.size - wind.size, mov.step) + wind.size - 1)) / 2) 
ind.fi <- (seq(1, crop.size - wind.size, mov.step) + wind.size - 1) + floor((crop.size - max(seq(1, crop.size - wind.size, mov.step) + wind.size - 1)) / 2)

#####
### CREATE AND WRITE RASTERS WITH LINEAR DETRENDING (NEW METRICS)
dir.mat <- matrix(data = NA, nrow = wind.size, ncol = wind.size)
for (i in 1:wind.size) {
  for (j in 1:wind.size) {
    if ((j - ((wind.size/2) + 1)) < 0) {
      dir.mat[i, j] <- atan((i - ((wind.size/2) + 1)) / (j - ((wind.size/2) + 1))) + pi
    } else {
      dir.mat[i, j] <- atan((i - ((wind.size/2) + 1)) / (j - ((wind.size/2) + 1)))
    }
  }
}
dir.mat <- dir.mat + pi/2

area <- readOGR(border.path)
area <- spTransform(area, crs(gt))
plot(area, add = T)
# if (prc*wind.size < 2) {part <- 2} else {part <- ceiling(prc*wind.size)}
part <- 0.05*wind.size

cl1 <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1)
registerDoParallel(cl1)
foreach (m = 1:tiles, .packages = c("tcltk", "raster", "spectral", "spatial", "minpack.lm")) %dopar% {
  if (!is.null(intersect(tt[m,], area))) { # 
    max.mag <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    max.imp <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    max.freq.lat <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    max.freq.long <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    pr.dir <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    pr.deg <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    wavelen <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    a_exp <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    l_exp <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    dev_exp <- matrix(NA, nrow = length(ind.st), ncol = length(ind.st))
    for (i in 1:length(ind.st)) {
      for (j in 1:length(ind.fi)) {
        ras.crop <- crop(gt.stack[[m]], extent(gt.stack[[m]], ind.st[i], ind.fi[i], ind.st[j], ind.fi[j]))
        # plot(ras.crop)
        data.sam <- as.matrix(ras.crop)
        if (sum(is.na(data.sam) | (data.sam == 0), max(table(data.sam))) < (wind.size^2 / 2) & !is.na(sum(data.sam))) {
          thr_dim <- sapply(as.data.frame(as.table(data.sam)), as.numeric)
          thr_dim <- thr_dim[!is.na(thr_dim[,3]),]
          topo_kr <- surf.ls(tr.deg, thr_dim[,1], thr_dim[,2], thr_dim[,3]) 
          zet <- matrix(predict.trls(topo_kr, x = rep(1:wind.size,wind.size), y = rep(1:wind.size, each = wind.size)), nrow = wind.size) 
          mat.detr <- data.sam - zet
          mat.detr[is.na(mat.detr)] <- 0
          four.im <- spec.fft(x = 1:wind.size, y = 1:wind.size, z = mat.detr)
          val <- unique(round(sort(abs(four.im$A), decreasing = T), di = 6))
          fi2 <- four.im
          fi2$A <- matrix(0, nrow = wind.size, ncol = wind.size)
          indc <- which(abs(four.im$A) > val[part+1])
          fi2$A[indc] <- four.im$A[indc]
          plt <- Re(spec.fft(fi2)$z)
          
          (max.imp[i,j] <- 1 - (sd(mat.detr - plt) / sd(mat.detr)))
          (max.mag[i,j] <- 4*max(abs(four.im$A)))
          
          dir <- dir.mat[indc]
          weit <- abs(four.im$A)[indc]
          xyt <- cbind(weit*sin(dir), weit*cos(dir))
          t <- prcomp(xyt)$rotation[1,1] / prcomp(xyt)$rotation[2,1]
          (pr.dir[i,j] <- 90 - 180*atan(t)/pi)
          (pr.deg[i,j] <- cor(xyt[,2] , xyt[,1])^2)
          
          if (!is.na(max(abs(four.im$A)))) {
            (max.freq.lat[i,j] <- abs((wind.size/2) - (which.max(abs(four.im$A)) %% wind.size) + 1))
            (max.freq.long[i,j] <- abs((wind.size/2) - floor(which.max(abs(four.im$A)) / wind.size)))
            (wavelen[i,j] <- (pix.size1*wind.size) / (max.freq.lat[i,j]^2 + max.freq.long[i,j]^2)^0.5)
          }
          
          if (NLR) {
            co <- mat.detr
            ro <- mat.detr
            for (ii in 1:wind.size) {
              co[,ii] <- ii
              ro[ii,] <- ii
            }
            cr <- (abs(co - (wind.size / 2 + 1))^2 + abs(ro - (wind.size / 2 + 1))^2)^0.5
            
            ud <- unique(sort(cr))
            ad <- ud
            for (ii in 1:length(ud)) {
              ad[ii] <- mean(abs(four.im$A)[which(cr == ud[ii])])
            }
            
            am <- abs(four.im$A)[match(unique(sort(cr)), cr)]
            dt <- cbind(ad, ud)[-1,] #[2:part]
            colnames(dt) <- c("ad", "ud")
            nlsm <- minpack.lm::nlsLM(formula = ad ~ a*exp(1)^(-l*ud), data = as.data.frame(dt), start = list(a = 1000, l = 1))
            (a_exp[i,j] <- coefficients(nlsm)[1])
            (l_exp[i,j] <- coefficients(nlsm)[2])
            (dev_exp[i,j] <- nlsm$m$deviance())
            
          }
          
        }
      }
      if(!exists("pb")) pb <- tcltk::tkProgressBar(title = paste0("Core # ",m), min=1, max=length(ind.st))
      tcltk::setTkProgressBar(pb, i)
      Sys.sleep(0.05)
      # print(i)
    }
    temp.ext <- tt.ext[[m]]
    temp.path <- file.path (folder.name,wind.size)
    max.mag.r <- raster(max.mag)
    max.imp.r <- raster(max.imp)
    max.freq.lat.r <- raster(max.freq.lat)
    max.freq.long.r <- raster(max.freq.long)
    pr.dir.r <- raster(pr.dir)
    pr.deg.r <- raster(pr.deg)
    wavelen.r <- raster(wavelen)
    a_exp.r <- raster(a_exp)
    l_exp.r <- raster(l_exp)
    dev_exp.r <- raster(dev_exp)
    extent(max.mag.r) <- extent(max.imp.r) <- extent(max.freq.lat.r) <- extent(max.freq.long.r) <- extent(pr.dir.r) <- extent(pr.deg.r) <-  extent(wavelen.r) <- extent(a_exp.r) <- extent(l_exp.r) <- extent(dev_exp.r) <- temp.ext
    writeRaster(max.mag.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_max.mag.tif"), format = "GTiff", overwrite=TRUE)
    writeRaster(max.imp.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_max.imp.tif"), format = "GTiff", overwrite=TRUE)
    # writeRaster(max.freq.lat.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_max.freq.lat.tif"), format = "GTiff", overwrite=TRUE)
    # writeRaster(max.freq.long.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_max.freq.long.tif"), format = "GTiff", overwrite=TRUE)
    writeRaster(pr.dir.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_pr.dir.tif"), format = "GTiff", overwrite=TRUE)
    writeRaster(pr.deg.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_pr.deg.tif"), format = "GTiff", overwrite=TRUE)
    writeRaster(wavelen.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_wavelen.tif"), format = "GTiff", overwrite=TRUE)
    writeRaster(a_exp.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_a_exp.tif"), format = "GTiff", overwrite=TRUE)
    writeRaster(l_exp.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_l_exp.tif"), format = "GTiff", overwrite=TRUE)
    writeRaster(dev_exp.r, filename = paste0(temp.path,"/", wind.size,"-",m,"_dev_exp.tif"), format = "GTiff", overwrite=TRUE)
  }
}
stopCluster(cl1)
print(paste0("The process is DONE with the window size of ", wind.size))
