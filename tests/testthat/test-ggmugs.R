test_that("ggmugs works", {
  plt = ggmugs(
    study_name = c("study1", "study2"),
    summary_stat = c("https://raw.githubusercontent.com/Broccolito/ggmugs_data/main/sumstat1.txt",
                     "https://raw.githubusercontent.com/Broccolito/ggmugs_data/main/sumstat2.txt"),
    p1 = 1e-4,
    p2 = 1e-6,
    p3 = 1e-8,
    color1 = "#FFFFE0",
    color2 = "#FFC300",
    color3 = "#FF5733"
  )
})
