#'Join tables
#'
#'@description Join tables resulting from other functions of the package \code{ntimes}.
#'
#'@param tab.x a data frame with attribute \code{ntimes}.
#'@param tab.y a data frame with attribute \code{ntimes}.
#'@param digits a numeric value specifying the number of digits for the p values.
#'@param alpha a numeric value specifying the significance level to indicate statistically significant multiple comparisons.
#'@param save a logical value indicating whether the output should be saved as a csv file.
#'@param file a character value indicating the name of output file in csv format to be saved.
#'@examples
#'library(dplyr)
#'library(magrittr)
#'
#'iris_nt <- iris %>% filter(Species != "versicolor") %>%
#'mutate(Species = droplevels(Species))
#'tab01 <- nt_describe(iris_nt, group = Species)
#'tab02 <- nt_compare_tg(iris_nt, group = Species)
#'tab <- nt_join_tables(tab01, tab02)
#'
#'@return A data frame containing \code{tab.x} and \code{tab.y}.
#'
#'@importFrom tidyr replace_na
#'@importFrom dplyr mutate select left_join
#'@importFrom stringr str_split
#'
#'@export
nt_join_tables <- function(tab.x, tab.y, digits = 3,
                         alpha = 0.05,
                         save = FALSE,
                         file = "table"){

  if ((attr(tab.x, "ntimes") == "descriptive" &
       (attr(tab.y, "ntimes") == "two_groups" |
        attr(tab.y, "ntimes") == "multiple_groups" |
        attr(tab.y, "ntimes") == "multiple_comparisons")) |
      attr(tab.y, "ntimes") == "descriptive" &
      (attr(tab.x, "ntimes") == "two_groups" |
       attr(tab.x, "ntimes") == "multiple_groups" |
       attr(tab.x, "ntimes") == "multiple_comparisons")){

    if (attr(tab.x, "ntimes") == "descriptive" &
        (attr(tab.y, "ntimes") == "two_groups" |
         attr(tab.y, "ntimes") == "multiple_groups" |
         attr(tab.y, "ntimes") == "multiple_comparisons")){
      descriptive.tab <- tab.x
      test.tab <- tab.y
    }

    if (attr(tab.y, "ntimes") == "descriptive" &
        (attr(tab.x, "ntimes") == "two_groups" |
         attr(tab.x, "ntimes") == "multiple_groups" |
         attr(tab.x, "ntimes") == "multiple_comparisons")){
      descriptive.tab <- tab.x
      test.tab <- tab.y
    }

    if (attr(test.tab, "ntimes") == "two_groups"){

      if (!any(names(test.tab) == "p value"))
        stop("'test.tab' does not have any column 'p value'.")

      test.tab <- test.tab %>% mutate(Variable = as.character(.data$Variable)) %>%
        select(.data$Variable, .data$`p value`)
      descriptive.tab <- descriptive.tab %>%
        mutate(Variable = as.character(.data$Variable))

      tab <- left_join(descriptive.tab, test.tab, by = "Variable")
      out <- tab %>% mutate(`p value` =
                              ifelse(.data$`p value` < 0.001, "< 0.001",
                                     round(.data$`p value`, digits)))  %>%
        replace_na(list(`p value` = ""))

      if (save)
        write.csv(out, file = paste0(file, ".csv"))
    } else if (attr(test.tab, "ntimes") == "multiple_groups"){

      test.tab <- test.tab$omnibus.test %>%
        mutate(Variable = as.character(.data$Variable)) %>%
        select(.data$Variable, .data$`p value`)
      descriptive.tab <- descriptive.tab %>%
        mutate(Variable = as.character(.data$Variable))

      tab <- left_join(descriptive.tab, test.tab, by = "Variable")
      out <- tab %>% mutate(`p value` =
                              ifelse(.data$`p value` < 0.001, "< 0.001",
                                     round(.data$`p value`, digits)))  %>%
        replace_na(list(`p value` = ""))

      if (save)
        write.csv(out, file = paste0(file, ".csv"))
    } else if (attr(test.tab, "ntimes") == "multiple_comparisons"){

      aux <- test.tab$omnibus.test %>%
        mutate(Variable = as.character(.data$Variable)) %>%
        select(.data$Variable, .data$`p value`)

      temp <- test.tab$mc.test %>% select(-.data$`95% CI`, -.data$Test, -.data$Group) %>%
        spread(key = .data$Hypothesis, value = "p value", drop = TRUE)

      for(i in 2:ncol(temp)){
        temp[, i] <- ifelse(temp[, i] < alpha, letters[i-1], "")
      }

      mc <- temp %>% unite(col = "Comparisons", -.data$Variable, sep = "")

      test.tab <- left_join(aux, mc, by = "Variable")
      tab <- left_join(descriptive.tab, test.tab, by = "Variable") %>%
        mutate(`p value` =
                 ifelse(.data$`p value` < 0.001, "< 0.001",
                        round(.data$`p value`, digits))) %>%
        replace_na(list(`p value` = "", Comparisons = ""))
      comparisons <- setNames(names(temp[-1]), letters[1:(ncol(temp)-1)])

      out <- list(tab = tab, comparisons = comparisons)

      if (save){
        write.csv(out, file = paste0(file, ".csv"))
        write.csv(comparisons, file = paste0(file, "_legend.csv"))
      }
    }
  }

  if (attr(tab.x, "ntimes") == "descriptive" & attr(tab.y, "ntimes") == "descriptive"){
    tab.x <- tab.x %>% mutate(id = 1:nrow(tab.x))
    tab.y <- tab.y %>% mutate(id = 1:nrow(tab.y))
    tab <- left_join(tab.x, tab.y, by = c("id", "Variable")) %>% select(-id)
    temp <- str_split(colnames(tab), ":")
    colnames(tab)[3:4] <- c(temp[[3]][2], temp[[4]][2])

    if (save){
      write.csv(out, file = paste0(file, ".csv"))
    }

    out <- tab
    attr(out, "ntimes") <- "descriptive"
  }

  return(out)
}

