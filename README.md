# functional_enrichment_app

Code for a toy shiny app that calculates functional enrichment, and allows interactive tweaking of gene list sizes, gene pathway list sizes e.t.c.

## Building a static (shinylive) version

This app can be exported to a static, serverless site using
[shinylive](https://posit-dev.github.io/r-shinylive/) — it compiles to
WebAssembly and runs entirely in the browser via webR, so the exported
`site/` directory can be hosted on any plain static file server (e.g.
GitHub Pages) with no Shiny server running anywhere.

1. Install R and the `shinylive` package:
   ```r
   install.packages("shinylive", repos = "https://cloud.r-project.org")
   ```

2. Make sure every package the app uses (`shiny`, `eulerr`, `scales`) is
   also installed locally — `shinylive::export()` determines what to
   bundle by inspecting locally installed packages, so they need to be
   present first:
   ```r
   install.packages(c("eulerr", "scales"), repos = "https://cloud.r-project.org")
   ```

3. Export the app to a `site/` directory:
   ```bash
   Rscript -e 'shinylive::export(appdir = ".", destdir = "site")'
   ```

4. Serve it locally to test, then open the printed URL in a browser:
   ```bash
   cd site && python3 -m http.server 8791
   ```
   First load is slow — webR and the bundled R packages need to download
   and initialize in the browser, which can take 10-30+ seconds.

This repo also has a GitHub Actions workflow
(`.github/workflows/deploy-shinylive.yml`) that runs this same export
automatically and publishes the result to GitHub Pages on every push to
`master`.

