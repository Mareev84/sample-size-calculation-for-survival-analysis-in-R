# =============================================================
# Расчёт размера выборки для анализа выживаемости (лог-ранк, Шёнфельд)
# Равномерный набор, административное окончание, учёт потерь к наблюдению
# =============================================================

## ------------------------------
## 1) ОСНОВНЫЕ ПАРАМЕТРЫ (изменяйте здесь)
## ------------------------------
alpha     <- 0.05        # уровень значимости (двусторонний)
power     <- 0.80        # мощность (обычно 0.8 или 0.9)
ratio <- 1   # соотношение n_t / n_c = 1:1 (n_t — терапия, n_c — контроль)

p_control <- 0.30        # 1-годичная доля событий в контроле  
RRR_set   <- c(0.10, 0.15, 0.20, 0.25)   # ожидаемое относительное снижение риска (RRR)

T_acc     <- 1.0         # длительность набора (лет), равномерный набор
T_follow  <- 0.0         # период дополнительного наблюдения после завершения набора (лет)

dropout_c_1y <- 0.10     # потери в контрольной группе, доля в год (0.10 = 10%/год)
dropout_t_1y <- 0.10     # потери в группе терапии, доля в год (0.10 = 10%/год)
# Для чувствительного анализа попробуйте 0.05 или 0.15; можно задать разные значения если ожидаются асимметричные потери

## ------------------------------
## 2) ПРЕДПОЛОЖЕНИЯ МОДЕЛИ
## ------------------------------
# - Набор пациентов равномерен во времени (U ~ Uniform(0, T_acc))
# - Пациенты, не потерянные к наблюдению, наблюдаются до административного окончания (T_end = T_acc + T_follow)
# - Потери к наблюдению (dropout) моделируются независимо от событий и приводят к раннему цензурированию
# - Модель выживаемости: экспоненциальная (S(t) = exp(-λ * t));
#   p_control (за 1 год) -> S_c(1) = 1 - p_control -> λ_c = -log(S_c(1))
# - Пропорциональные риски: λ_t = HR * λ_c, где HR = 1 - RRR
# - Ожидаемая доля наблюдаемых событий (с учётом равномерного набора и потерь):
#   P(event) = (λ / (λ+δ)) * [1 - exp(-(λ+δ)*T_follow) * (1 - exp(-(λ+δ)*T_acc)) / ((λ+δ)*T_acc)]
# - Требуемое число событий по Шёнфельду:
#   D = ((z_{1-α/2} + z_{power})^2) / (log(HR))^2
# - Перевод числа событий в общий размер выборки:
#   n_total = D * (1 + r) / (p_obs_c + r * p_obs_t)

## ------------------------------
## 3) ФУНКЦИИ
## ------------------------------
required_events_schoenfeld <- function(HR, alpha = 0.05, power = 0.80) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)
  ((z_alpha + z_beta)^2) / (log(HR)^2)
}

lambda_control_from_p1y <- function(p_control_1y) {
  S1 <- 1 - p_control_1y
  -log(S1)
}

delta_from_dropout <- function(dropout_1y) {
  if (dropout_1y <= 0) return(0)
  -log(1 - dropout_1y)
}

# Вероятность события при равномерном наборе, админ. окончании и потерях
p_event_uniform_with_dropout <- function(lambda, delta, T_acc, T_follow) {
  lam <- lambda + delta
  if (T_acc == 0) {
    return((lambda/lam) * (1 - exp(-lam * T_follow)))
  }
  (lambda/lam) * (1 - exp(-lam * T_follow) * (1 - exp(-lam * T_acc)) / (lam * T_acc))
}

# Главная функция: расчёт общего n с учётом потерь в обеих группах
required_n_with_uniform_accrual_dropout <- function(p_control_1y, RRR,
                                                    T_acc, T_follow,
                                                    dropout_c_1y = 0.10,
                                                    dropout_t_1y = 0.10,
                                                    alpha = 0.05, power = 0.80, r = 1) {
  HR <- 1 - RRR
  D  <- required_events_schoenfeld(HR, alpha, power)
  
  lambda_c <- lambda_control_from_p1y(p_control_1y)
  lambda_t <- lambda_c * HR
  
  delta_c <- delta_from_dropout(dropout_c_1y)
  delta_t <- delta_from_dropout(dropout_t_1y)
  
  p_obs_c <- p_event_uniform_with_dropout(lambda_c, delta_c, T_acc, T_follow)
  p_obs_t <- p_event_uniform_with_dropout(lambda_t, delta_t, T_acc, T_follow)
  
  n_total <- D * (1 + r) / (p_obs_c + r * p_obs_t)
  ceiling(n_total)
}

## ------------------------------
## 4) РАСЧЁТ ТАБЛИЦЫ
## ------------------------------
n_vec <- sapply(RRR_set, function(rrr)
  required_n_with_uniform_accrual_dropout(p_control, rrr, T_acc, T_follow,
                                          dropout_c_1y = dropout_c_1y,
                                          dropout_t_1y = dropout_t_1y,
                                          alpha = alpha, power = power, r = ratio)
)

tbl <- data.frame(
  `Снижение риска (RRR, %)` = round(RRR_set * 100),
  `Общее n`                 = n_vec,
  `≈ n на группу`           = ceiling(n_vec / (1 + ratio))
)
print(tbl, row.names = FALSE)

## ------------------------------
## 5) ГРАФИК (с подписями n)
## ------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

df_plot <- data.frame(RRR_pct = round(RRR_set * 100), n = n_vec)

title_ru <- sprintf(
  "Размер выборки при p_control=%d%%; равномерный набор %.1f г., доп. наблюдение %.1f г.; потери=%.0f%%/год\n(α=%.2f, мощность=%.0f%%, 1:%d; экспоненциальная модель)",
  round(p_control*100), T_acc, T_follow, dropout_c_1y*100, alpha, power*100, ratio+1
)

p <- ggplot(df_plot, aes(x = RRR_pct, y = n)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = n), df_plot, linetype = "dotted", linewidth = 0.5) +
  geom_text(aes(label = n), vjust = -0.6, size = 3.8) +
  scale_x_continuous(breaks = df_plot$RRR_pct) +
  labs(
    title = title_ru,
    x = "Относительное снижение риска (RRR), %",
    y = "Требуемый общий размер выборки (n)"
  ) +
  theme_bw()

print(p)

## ------------------------------
## (опционально) Сохранение результатов
## ------------------------------
# write.csv(tbl, "sample_size_with_dropout_table.csv", row.names = FALSE)
# ggsave("sample_size_with_dropout_plot.png", p, width = 8, height = 5, dpi = 150)


