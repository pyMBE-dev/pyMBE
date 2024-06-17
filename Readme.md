# Corporate design

The main logo is [logo.svg](logo.svg).

All derived logos can be regenerated with `make all`.
To create a new logo derived from the main one, add an extra target in the Makefile.
Use [Inkscape Actions](https://gist.github.com/jasonm23/c323bc20380d5569de68de23ee8d2ffa)
to transform individual elements of the SVG file.

## Continuous deployment

The generated files are automatically deployed to branch `logos`. This is achieved
using [GitHub Pages action](https://github.com/marketplace/actions/github-pages-action)
with a deployment token[^peaceiris-actions-gh-pages-faq-ssh].

To regenerate the deploy key and deployment secret, follow these steps:

1. create a new passwordless public/private key pair with this command[^github-docs-ssh]:
   ```sh
   ssh-keygen -t ed25519 -f gh-pages -N "" -C "41898282+github-actions[bot]@users.noreply.github.com"
   cat ./gh-pages.pub
   cat ./gh-pages
   ```
2. create a new repository deploy key to store the public key[^github-docs-deploy-keys]
   * got to the repository settings
   * click "Deploy keys"
   * click "Add deploy key"
   * write `deploy-logos` in the "Title" field
   * copy-paste the contents of `./gh-pages.pub` in the "Key" field
   * check the "Allow write access" box
   * click "Add key".
3. create a new repository environment to store the private key[^github-docs-env-secrets]
   * go to the repository settings
   * click "Environments"
   * click "New environments"
   * enter name `deploy-logos`
   * click "Configure environment"
   * in "Deployment branches and tags", click "No restriction" and "Selected branches and tags"
   * click "Add deployment branch or tag rule"
   * enter name pattern `corporate-design`
   * click "Add rule"
   * click "Add environment secret"
   * write `ACTIONS_DEPLOY_KEY` in the "Name" field
   * copy-paste the contents of `./gh-pages` in the "Value" field
   * click "Add secret"

[^peaceiris-actions-gh-pages-faq-ssh]: FAQ section [Create SSH Deploy Key](https://github.com/peaceiris/actions-gh-pages/tree/v4.0.0?tab=readme-ov-file#tips-and-faq) explains the steps to register a SSH public/private key pair.
[^github-docs-ssh]: GitHub documentation: [Generating a new SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key) shows the exact command line options, in particular which algorithm to choose.
[^github-docs-deploy-keys]: GitHub documentation: [Set up deploy keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/managing-deploy-keys#set-up-deploy-keys)
[^github-docs-env-secrets]: GitHub documentation: [Creating secrets for an environment](https://docs.github.com/en/actions/security-guides/using-secrets-in-github-actions#creating-secrets-for-an-environment)
