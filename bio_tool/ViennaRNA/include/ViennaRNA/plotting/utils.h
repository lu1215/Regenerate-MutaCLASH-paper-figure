#ifndef VIENNA_RNA_PACKAGE_PLOT_UTILS_H
#define VIENNA_RNA_PACKAGE_PLOT_UTILS_H

/**
 *  @file     ViennaRNA/plotting/utils.h
 *  @ingroup  utils, plotting_utils
 *  @brief    Various utilities to assist in plotting secondary structures and consensus structures
 */

/**
 *  @addtogroup  plotting_utils
 *  @{
 *  @brief  Functions for Creating Secondary Structure Plots, Dot-Plots, and More
 *  @}
 */


/**
 *  @addtogroup  annotation_utils
 *  @{
 *  @brief  Functions to generate annotations for Secondary Structure Plots, Dot-Plots, and Others
 */

/**
 *  @brief  Produce covariance annotation for an alignment given a secondary structure
 *
 */
char **
vrna_annotate_covar_struct(const char **alignment,
                           const char *structure,
                           vrna_md_t  *md);


/**
 *  @brief  Produce covariance annotation for an alignment given a set of base pairs
 *
 */
vrna_cpair_t *
vrna_annotate_covar_pairs(const char  **alignment,
                          vrna_ep_t   *pl,
                          vrna_ep_t   *mfel,
                          double      threshold,
                          vrna_md_t   *md);


/**
 * @}
 */

#endif
