/*
Language: YAML
Author: Jan Rezac
*/

hljs.LANGUAGES.yaml =
{
  case_insensitive: true,
  defaultMode: {
    contains: ['comment',  'number', 'string', 'label'],
    illegal: ''
  },

  modes: [
    hljs.C_BLOCK_COMMENT_MODE,
    {
      className: 'comment',
      begin: '#',  end: '$'
    },
    hljs.BACKSLASH_ESCAPE,
    hljs.QUOTE_STRING_MODE,
    {
      className: 'string',
      begin: '\'', end: '[^\\\\]\'',
      illegal: '[^\\\\][^\']'
    },
    hljs.C_NUMBER_MODE,
    {
      className: 'label',
      begin: '^\\s*[A-Za-z0-9_.$]+:',  end: '^'
    },
    {
      className: 'label',
      begin: '^\\s*-\\s',  end: '^'
    },
  ]
};

