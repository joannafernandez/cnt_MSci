# Scientific Graphing: Best Practices for Clear Communication

## Importnat resources
https://urbaninstitute.github.io/graphics-styleguide/
https://practicaldatascience.org/notebooks/class_5/week_1/2.1.1_the_matplotlib_model.html

https://www.nature.com/articles/s41467-020-19160-7
https://rupress.org/jcb/article/219/6/e202001064/151717/SuperPlots-Communicating-reproducibility-and

colours:

https://carto.com/carto-colors/

https://cran.r-project.org/web/packages/khroma/vignettes/tol.html#sec:light

https://www.kennethmoreland.com/color-maps/ColorMapsExpanded.pdf

## 🎯 Core Principle
**A figure is not decoration — it is an argument.**  
Every graph should communicate one clear message.

- Define the *take-home message* before plotting
- Design every element to support that message
- Avoid unnecessary complexity or visual noise

---

## 1. Show the Data, Not Just the Summary

### Key idea (from SuperPlots paper)
- Traditional bar plots hide variability and reproducibility
- Always aim to show **raw data + summary statistics**

### Good practice
- Plot individual data points (e.g., scatter, jitter, violin)
- Overlay summary (mean ± SD/SEM)
- Distinguish:
  - **technical replicates**
  - **biological replicates**

### Recommended formats
- Dot plots
- Box/violin + points
- “SuperPlots” (group-level + replicate-level data)

> Showing underlying data improves transparency and reproducibility.

---

## 2. Be Honest and Avoid Misleading Visuals

- Do not manipulate axes to exaggerate effects  
- Use **zero baseline** where appropriate (especially bar plots)
- Keep aspect ratios sensible

> Misleading graphs damage scientific credibility

---

## 3. Choose the Right Plot Type

### Prefer simple, interpretable plots
- Line plots → trends
- Scatter plots → relationships
- Bar plots → simple comparisons (with caution)
- Box/violin → distributions

Avoid:
- Overly complex or novel plots unless necessary  
- Pie charts for precise comparisons

> Simple plots are often more effective than complex ones :contentReference[oaicite:0]{index=0}

---

## 4. Label Everything Clearly

- Axes must include:
  - variable name
  - units
- Use human-readable labels (not variable names)
- Include:
  - legend
  - descriptive title or caption

> No ambiguity — the figure should stand alone

---

## 5. Use Colour Thoughtfully

### From Urban Institute style guide
- Use colour to **encode meaning**, not decorate :contentReference[oaicite:1]{index=1}
- Keep palettes:
  - consistent
  - minimal

### Good practice
- Use colourblind-friendly palettes
- Avoid red–green contrasts
- Do not rely on colour alone (add shapes/labels)

---

## 6. Reduce Visual Clutter (Maximise Data-to-Ink Ratio)

- Remove:
  - unnecessary gridlines
  - background shading
  - excessive borders
- Keep only elements that add information

> “Keep it simple” improves readability :contentReference[oaicite:2]{index=2}

---

## 7. Make Figures Accessible

- Use readable font sizes
- Ensure sufficient contrast
- Avoid tiny annotations

From Urban guidelines:
- Typography should create **visual hierarchy** (title > axes > labels) :contentReference[oaicite:3]{index=3}

---

## 8. Represent Variability and Uncertainty Properly

- Always show:
  - SD, SEM, or confidence intervals
- Clearly state what error bars represent
- Avoid hiding variability with only means

---

## 9. Ensure Reproducibility is Visible

### From SuperPlots concept
- Highlight independent replicates explicitly
- Avoid pooling everything into a single distribution

Example:
- Colour points by replicate
- Show replicate means

---

## 10. Maintain Consistency Across Figures

- Same:
  - axis scales
  - colours
  - symbols
- Enables easy comparison between panels

---

## 11. Know Your Audience

- Tailor complexity to reader expertise
- Use annotations if needed to guide interpretation

> A plot should match the audience’s ability to interpret it :contentReference[oaicite:4]{index=4}

---

## 12. Provide Context in the Caption

A good caption should include:
- What is being shown
- Experimental conditions
- Sample size (n)
- Statistical test used

---

## 13. Test Your Figure

Before finalising:
- Can someone understand it in <10 seconds?
- Is the main message obvious?
- Is anything misleading?

---

## ✅ Quick Checklist

- [ ] Clear message
- [ ] Raw data shown
- [ ] Appropriate plot type
- [ ] Axes labelled with units
- [ ] Colour used meaningfully
- [ ] Minimal clutter
- [ ] Variability shown
- [ ] Replicates visible
- [ ] Consistent formatting
- [ ] Caption explains everything

---

## 📚 Key References

- Nature Communications (2020): Best practices in data visualization  
- SuperPlots: *Communicating reproducibility and variability in cell biology*  
- Urban Institute Data Visualization Style Guide :contentReference[oaicite:5]{index=5}  
- General plotting principles :contentReference[oaicite:6]{index=6}  

---

## 💡 Final Thought

> “Clarity beats cleverness.”  
If a figure needs explaining, it’s not finished.
