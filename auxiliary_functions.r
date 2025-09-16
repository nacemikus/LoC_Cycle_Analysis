
logit <- function(x) log(x/(1-x))

inv_logit <- function(x) 1/(1+exp(-x))

pw = function(x,delta) x^delta/(x^delta + (1-x)^delta) # prob weighting


# inv_logit2 <- function(x,a=1) 1/(1+exp(-a*x))

phi_approx <- function(x) inv_logit(0.07056*x^3 + 1.5976*x)

# get CrIs for text
sf <- function(x, pub = 0, prob_vect = c(0.5,0.025,0.975), dec_no = 3) { 
  y <- quantile(x, probs=prob_vect)%>% round(dec_no) 
  y[[4]] <- mean(x<0) %>% round(3)
  names(y)[4] <- "p"
  p_val = y[[4]];
  if (y[[4]] > 0.5) y[[4]] = 1 - y[[4]]
  
  if (y[[1]] < 0) {
    y_text = paste("b = ", y[[1]], ", 95% CrI [",y[[2]], ", ",y[[3]], "], P(b>0) = ",y[[4]], sep ="")
  } else   y_text = paste("b = ", y[[1]], ", 95% CrI [",y[[2]], ", ",y[[3]], "], P(b<0) = ",y[[4]], sep ="")
  
  if (pub == 0) {
    return(y)
  } else {
    
    return(y_text)  
  } }

# remove outliers
wo <- function(x, d = 3, remove.na = FALSE) {
  if (remove.na) {
    x <- x[abs(x - mean(x, na.rm =TRUE)) <d*sd(x, na.rm =TRUE)]
  } else  {
    x[abs(x - mean(x, na.rm =TRUE)) > d*sd(x, na.rm =TRUE)] <- NA
    
  }
  return(x)}
