# Analytical figures for genes8PT. Base R only; no stochastic simulations.
# Rscript genes8PT_analytical_figures.R [output_directory]
# Figures S3/S4 reproduce the main figures in centimorgans (cM).
# Main text: recombination only. Mutation appears only in Figures S1/S2.
# Frequency-dependent persistence uses a continuous two-class model in which
# breakup slows allele spread. It is not a general multilocus simulation.
args <- commandArgs(trailingOnly=TRUE)
out_dir <- if(length(args)) args[1] else "genes8PT_analytical_figures"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
Ne_values <- c(10,50,200,1000)
s_min <- 0.001
s_max <- 0.25
n_curve <- 2000L
rho <- 1e-8 # Morgan/base pair: illustrative 1 cM/Mb, not a universal map
break_mode <- "endpoint" # alternatively "any_crossover"
main1_max_bp <- 30e6
main2_max_bp <- 3e6
mu <- 1.2e-8 # disruptive mutation probability per required base/generation
make_supplementary <- TRUE
stopifnot(rho>0,mu>=0,mu<1,all(Ne_values>1),s_min>0,s_max>s_min,
          break_mode %in% c("endpoint","any_crossover"))

# B is inclusive physical span; B-1 is the separation between endpoints.
map_length <- function(B,rho_value=rho) rho_value*(B-1)
recombination_loss <- function(B) {
  stopifnot(all(B>=1))
  L <- map_length(B)
  if(break_mode=="endpoint") -expm1(-2*L)/2 else -expm1(-L)
}
inverse_recombination <- function(c) {
  limit <- if(break_mode=="endpoint") 0.5 else 1
  answer <- rep(Inf,length(c))
  answer[c<=0] <- 1
  valid <- c>0 & c<limit
  L <- if(break_mode=="endpoint") -log1p(-2*c[valid])/2 else -log1p(-c[valid])
  answer[valid] <- 1+L/rho
  answer
}
mutation_loss <- function(B) -expm1(B*log1p(-mu))
total_loss <- function(B,mutation=FALSE) {
  c <- recombination_loss(B)
  if(mutation) 1-(1-c)*(1-mutation_loss(B)) else c
}
effective_selection <- function(s,B,mutation=FALSE) (1+s)*(1-total_loss(B,mutation))-1
half_sweep <- function(s,Ne,B,mutation=FALSE) {
  lc <- -log1p(-recombination_loss(B))
  lu <- if(mutation) -B*log1p(-mu) else rep(0,length(B))
  b <- log1p(s)-lc
  a <- b-lu
  x0 <- 1/(2*Ne)
  valid <- a>0 & b>0 & a-b/2>0
  T <- A <- rep(Inf,length(B))
  T[valid] <- (log((0.5/x0)*((a[valid]-b[valid]*x0)/(a[valid]-b[valid]/2))))/a[valid]
  A[valid] <- (log(Ne)+lu[valid]*T[valid])/b[valid]
  list(T=T,A=A,reachable=valid)
}
log_retention <- function(s,Ne,B,mutation=FALSE) {
  times <- half_sweep(s,Ne,B,mutation)
  result <- rep(-Inf,length(B))
  valid <- times$reachable
  result[valid] <- log1p(-recombination_loss(B[valid]))*times$A[valid]
  if(mutation) result[valid] <- result[valid]+B[valid]*log1p(-mu)*times$T[valid]
  result
}

# Recombination-only boundaries are explicit inverse formulas.
selection_bound <- function(s,Ne,mutation=FALSE,max_bp=main1_max_bp) {
  cm <- (s-1/Ne)/(1+s)
  if(!mutation) return(ifelse(cm>0,inverse_recombination(cm),0))
  vapply(s,function(si) {
    if(effective_selection(si,1,TRUE)<=1/Ne) return(0)
    f <- function(B) effective_selection(si,B,TRUE)-1/Ne
    if(f(max_bp)>=0) return(max_bp)
    uniroot(f,c(1,max_bp),tol=1e-5)$root
  },numeric(1))
}
retention_bound <- function(s,Ne,p,mutation=FALSE,max_bp=main2_max_bp) {
  if(!mutation) {
    K <- -log(p)/log(Ne)
    return(inverse_recombination(-expm1(-log1p(s)*K/(1+K))))
  }
  vapply(s,function(si) {
    f <- function(B) {
      lp <- log_retention(si,Ne,B,TRUE)
      if(is.finite(lp)) lp-log(p) else -1e6
    }
    if(f(1)<0) return(0)
    if(f(max_bp)>=0) return(max_bp)
    uniroot(f,c(1,max_bp),tol=1e-5)$root
  },numeric(1))
}
# For a mutation extension, these scalar roots locate boundaries of explicit
# analytical functions; they do not simulate allele trajectories.

axis_value <- function(B,unit) if(unit=="cM") 100*map_length(pmax(1,B)) else B/1e6

left_region <- function(bound,s,max_bp,col,unit) {
  polygon(c(0,axis_value(pmax(0,pmin(max_bp,bound)),unit),0),
          c(s_min,s,s_max),col=col,border=NA)
}
right_region <- function(bound,s,max_bp,col,unit) {
  polygon(axis_value(c(max_bp,pmax(0,pmin(max_bp,bound)),max_bp),unit),
          c(s_min,s,s_max),col=col,border=NA)
}
draw <- function(number,mutation=FALSE,unit="Mb") {
  s <- seq(s_min,s_max,length.out=n_curve)
  max_bp <- if(number==1) main1_max_bp else main2_max_bp
  # Keep the genetic-axis ranges fixed when physical map calibration changes.
  if(unit=="cM") max_bp <- 1+(if(number==1) 0.30 else 0.03)/rho
  cols <- c("#2C7FB8","#F0E17A","#D9793D","#E0E0E0")
  layout(matrix(c(1,2,3,4,5,5),nrow=3,byrow=TRUE),heights=c(1,1,0.48))
  par(mar=c(4.8,5.5,2.1,1.0),oma=c(0,0,3,0),
      mgp=c(3.6,0.9,0),las=1,family="sans",cex.axis=1.25,cex.lab=1.25,cex.main=1.35)
  for(Ne in Ne_values) {
    plot(NA,xlim=c(0,axis_value(max_bp,unit)),ylim=c(s_min,s_max),xaxs="i",yaxs="i",
         xlab=if(unit=="cM") "Genetic map distance (cM)" else "B, physical span (Mb)",ylab="s, selective advantage",
         main=bquote(N[e]==.(Ne)))
    rect(0,s_min,axis_value(max_bp,unit),s_max,col=cols[3],border=NA)
    if(number==1) {
      # At a fixed s, the lower negative boundary has larger admissible B.
      if(!mutation) {
        negative <- inverse_recombination((s+1/Ne)/(1+s))
      } else {
        negative <- vapply(s,function(si) {
          f <- function(B) effective_selection(si,B,TRUE)+1/Ne
          if(f(max_bp)>=0) return(max_bp)
          uniroot(f,c(1,max_bp),tol=1e-5)$root
        },numeric(1))
      }
      positive <- selection_bound(s,Ne,mutation,max_bp)
      left_region(negative,s,max_bp,cols[2],unit)
      left_region(positive,s,max_bp,cols[1],unit)
    } else {
      b50 <- retention_bound(s,Ne,0.50,mutation,max_bp)
      b95 <- retention_bound(s,Ne,0.95,mutation,max_bp)
      left_region(b50,s,max_bp,cols[2],unit)
      left_region(b95,s,max_bp,cols[1],unit)
      selectable <- selection_bound(s,Ne,mutation,max_bp)
      if(mutation) {
        reachable <- vapply(s,function(si) {
          f <- function(B) log1p(si)+log1p(-recombination_loss(B))+2*B*log1p(-mu)
          if(f(1)<=0) return(0)
          if(f(max_bp)>=0) return(max_bp)
          uniroot(f,c(1,max_bp),tol=1e-5)$root
        },numeric(1))
        selectable <- pmin(selectable,reachable)
      }
      # The 95% and 50% boundaries are smooth analytical curves, same style.
      for(bound in list(b50,b95)) {
        visible <- bound>=1 & bound<=max_bp & bound<selectable
        segments <- rle(visible)
        end <- cumsum(segments$lengths)
        start <- end-segments$lengths+1
        for(i in which(segments$values)) {
          j <- start[i]:end[i]
          lines(axis_value(bound[j],unit),s[j],lwd=0.6)
        }
      }
      right_region(selectable,s,max_bp,cols[4],unit)
    }
    box()
  }
  name <- if(unit=="cM") paste0("Figure S",number+2) else
    paste0("Figure ",if(mutation) "S" else "",number)
  title <- if(number==1) if(unit=="cM") "Initial selection and genetic map distance" else "Initial selection and physical span" else
    "Analytical retention during the initial half sweep"
  mtext(paste(name,title),outer=TRUE,side=3,line=1.5,font=2,cex=1.05)
  calibration <- paste0("Illustrative map: ",format(rho*1e8)," cM/Mb; ",break_mode,
                        if(mutation) "; mutation extension" else "; recombination only")
  if(unit=="cM") calibration <- "Haldane endpoint map; recombination only; 100 cM = 1 Morgan"
  mtext(calibration,outer=TRUE,side=3,line=0.2,cex=0.85)
  par(mar=c(0,0,0,0))
  plot.new()
  labels <- if(number==1) c("Positive net selection","Drift zone","Negative net selection") else
    c("P > 0.95","0.50 <= P <= 0.95","P < 0.50",
      if(mutation) "Initial criterion / half sweep not met" else "Initial criterion not met")
  legend("center",legend=labels,fill=cols[seq_along(labels)],ncol=1,bty="n",cex=1.9,x.intersp=0.7,y.intersp=1.05)
}
save_fig <- function(number,mutation=FALSE,unit="Mb") {
  stem <- file.path(out_dir,paste0(if(unit=="cM") paste0("FigureS",number+2) else
    paste0("Figure",if(mutation) "S" else "",number),"_analytical"))
  pdf(paste0(stem,".pdf"),width=9,height=8.5,useDingbats=FALSE)
  tryCatch(draw(number,mutation,unit),finally=dev.off())
  png(paste0(stem,".png"),width=2700,height=2550,res=300)
  tryCatch(draw(number,mutation,unit),finally=dev.off())
}

# Verify closed-form thresholds against independent defining expressions.
for(Ne in Ne_values) {
  s <- c(0.11,0.17,0.25)
  B <- selection_bound(s,Ne)
  stopifnot(max(abs(effective_selection(s,B)-1/Ne))<1e-12)
  for(p in c(0.50,0.95)) {
    B <- retention_bound(s,Ne,p)
    stopifnot(max(abs(exp(log_retention(s,Ne,B))-p))<1e-12)
  }
  # An independent quadrature verifies the analytical frequency exposure.
  si <- 0.17
  x0 <- 1/(2*Ne)
  Btest <- 125001
  gamma <- log1p(si)+log1p(-recombination_loss(Btest))
  times <- half_sweep(si,Ne,Btest)
  numerical <- integrate(function(t) 1-plogis(qlogis(x0)+gamma*t),0,times$T)$value
  stopifnot(abs(numerical-times$A)<1e-8)
}
human <- expand.grid(s=c(0.001,0.01,0.05),Ne=c(1000,10000,100000))
human$initial_bound_bp <- mapply(function(s,Ne) selection_bound(s,Ne),human$s,human$Ne)
human$B95_bp <- mapply(function(s,Ne) retention_bound(s,Ne,0.95),human$s,human$Ne)
human$B50_bp <- mapply(function(s,Ne) retention_bound(s,Ne,0.50),human$s,human$Ne)
save_fig(1)
save_fig(2)
save_fig(1,unit="cM")
save_fig(2,unit="cM")
if(make_supplementary) {
  save_fig(1,TRUE)
  save_fig(2,TRUE)
}
write.csv(human,file.path(out_dir,"human_recombination_only_bounds.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out_dir,"R_session_info.txt"))
cat("Analytical figures created in",normalizePath(out_dir),"\n")